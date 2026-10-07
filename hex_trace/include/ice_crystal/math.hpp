#pragma once

#include <cmath>
#include <cfloat>
#include <cstdint>

#if defined(__CUDACC__)
#define ICE_HD __host__ __device__
#else
#define ICE_HD
#endif

namespace iceCrystal
{
inline constexpr double pi = 3.141592653589793238462643383279502884;
inline constexpr double sqrt3 = 1.732050807568877293527446341505872367;

ICE_HD inline bool finite_value(double x) noexcept { return x >= -DBL_MAX && x <= DBL_MAX; }
ICE_HD inline double max_value(double a, double b) noexcept { return a > b ? a : b; }

// Independent FP64 types: the existing rain-drop field and interface remain FP32
// and retain their electric-field (not ray-flux) contract.
struct Vec3
{
    double x = 0, y = 0, z = 0;
    ICE_HD Vec3 operator+(Vec3 b) const noexcept { return {x+b.x, y+b.y, z+b.z}; }
    ICE_HD Vec3 operator-(Vec3 b) const noexcept { return {x-b.x, y-b.y, z-b.z}; }
    ICE_HD Vec3 operator-() const noexcept { return {-x,-y,-z}; }
    ICE_HD Vec3 operator*(double b) const noexcept { return {x*b,y*b,z*b}; }
    ICE_HD Vec3 operator/(double b) const noexcept { return {x/b,y/b,z/b}; }
    ICE_HD double dot(Vec3 b) const noexcept { return x*b.x+y*b.y+z*b.z; }
    ICE_HD Vec3 cross(Vec3 b) const noexcept
    { return {y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x}; }
    ICE_HD bool is_finite() const noexcept { return finite_value(x)&&finite_value(y)&&finite_value(z); }
    ICE_HD double max_abs() const noexcept
    { return max_value(::fabs(x),max_value(::fabs(y),::fabs(z))); }
};
ICE_HD inline bool normalize(Vec3 a, Vec3& result) noexcept
{
    const double s=a.max_abs();
    if(!(s>0) || !finite_value(s)) return false;
    a=a/s;
    const double n=::sqrt(a.dot(a));
    if(!(n>0) || !finite_value(n)) return false;
    result=a/n;
    return result.is_finite();
}
ICE_HD inline bool is_unit(Vec3 a) noexcept
{ return a.is_finite() && ::fabs(a.dot(a)-1.0)<=128*DBL_EPSILON; }
struct Frame
{
    Vec3 e0, e1;
    ICE_HD bool valid(Vec3 k) const noexcept
    {
        return is_unit(k)&&is_unit(e0)&&is_unit(e1) &&
            ::fabs(k.dot(e0))<=128*DBL_EPSILON &&
            ::fabs(k.dot(e1))<=128*DBL_EPSILON &&
            (k.cross(e0)-e1).max_abs()<=256*DBL_EPSILON;
    }
};
ICE_HD inline bool make_frame(Vec3 k, Frame& f) noexcept
{
    if(!is_unit(k)) return false;
    const double ax=::fabs(k.x),ay=::fabs(k.y),az=::fabs(k.z);
    const Vec3 a=(ax<=ay && ax<=az)?Vec3{1,0,0}:(ay<=az?Vec3{0,1,0}:Vec3{0,0,1});
    if(!normalize(a-k*k.dot(a),f.e0)) return false;
    if(!normalize(k.cross(f.e0),f.e1)) return false;
    return f.valid(k);
}
struct Complex
{
    double re=0, im=0;
    ICE_HD Complex operator+(Complex b) const noexcept { return {re+b.re,im+b.im}; }
    ICE_HD Complex operator*(double b) const noexcept { return {re*b,im*b}; }
    ICE_HD Complex operator*(Complex b) const noexcept
    { return {re*b.re-im*b.im,re*b.im+im*b.re}; }
    ICE_HD double norm2() const noexcept { return re*re+im*im; }
};
struct Field { Complex x,y; };
struct JonesFlux
{
    // Columns: responses to the same two orthogonal unit incident fields.
    // Rows: components in the current transverse frame. NOT a single Jones
    // vector purporting to represent unpolarized light.
    Field column[2]{};
    ICE_HD static JonesFlux identity() noexcept
    {
        JonesFlux j{}; j.column[0].x.re=1; j.column[1].y.re=1; return j;
    }
    ICE_HD double power() const noexcept
    {
        return 0.5*(column[0].x.norm2()+column[0].y.norm2()+
                    column[1].x.norm2()+column[1].y.norm2());
    }
    ICE_HD bool is_finite() const noexcept
    {
        for(unsigned i=0;i<2;++i)
            if(!finite_value(column[i].x.re)||!finite_value(column[i].x.im)||
               !finite_value(column[i].y.re)||!finite_value(column[i].y.im)) return false;
        return true;
    }
    ICE_HD bool nonzero() const noexcept
    {
        for(unsigned i=0;i<2;++i)
            if(column[i].x.re!=0 || column[i].x.im!=0 ||
               column[i].y.re!=0 || column[i].y.im!=0) return true;
        return false;
    }
};
ICE_HD inline JonesFlux change_frame(const JonesFlux& j, Frame old, Frame next) noexcept
{
    JonesFlux r{};
    for(unsigned c=0;c<2;++c)
    {
        r.column[c].x=j.column[c].x*old.e0.dot(next.e0)+j.column[c].y*old.e1.dot(next.e0);
        r.column[c].y=j.column[c].x*old.e0.dot(next.e1)+j.column[c].y*old.e1.dot(next.e1);
    }
    return r;
}
ICE_HD inline JonesFlux scale_rows(JonesFlux j, Complex s, Complex p) noexcept
{
    for(unsigned c=0;c<2;++c) { j.column[c].x=j.column[c].x*s; j.column[c].y=j.column[c].y*p; }
    return j;
}
struct Sum
{
    double value=0, correction=0;
    ICE_HD void add(double x) noexcept
    {
        const double y=x-correction,t=value+y;
        correction=(t-value)-y; value=t;
    }
};
// Stateless, reproducible counter scrambling. Midpoints of 2^52 cells in (0,1).
// This is an input quadrature/MC source, not an optical approximation.
ICE_HD inline std::uint64_t mix64(std::uint64_t x) noexcept
{
    x=(x^(x>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    x=(x^(x>>27))*UINT64_C(0x94d049bb133111eb);
    return x^(x>>31);
}
ICE_HD inline double sample_uniform(std::uint64_t seed,std::uint64_t id,unsigned dimension) noexcept
{
    const auto x=mix64(seed+UINT64_C(0x9e3779b97f4a7c15)*(3*id+dimension+1));
    return (double(x>>12)+0.5)*0x1p-52;
}
} // namespace iceCrystal
