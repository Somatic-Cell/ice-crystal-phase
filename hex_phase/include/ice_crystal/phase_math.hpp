#pragma once
#include <ice_crystal/hex_trace_core.hpp>

namespace iceCrystal
{
struct Rotation
{
    Vec3 x{1,0,0},y{0,1,0},z{0,0,1}; // columns: body -> ensemble reference frame
    ICE_HD Vec3 apply(Vec3 p) const noexcept {return x*p.x+y*p.y+z*p.z;}
    ICE_HD Vec3 inverse(Vec3 p) const noexcept {return {x.dot(p),y.dot(p),z.dot(p)};}
    ICE_HD bool valid() const noexcept
    {
        return is_unit(x)&&is_unit(y)&&is_unit(z)&&
            ::fabs(x.dot(y))<=128*DBL_EPSILON&&
            (x.cross(y)-z).max_abs()<=256*DBL_EPSILON;
    }
};
// Scalar order w,x,y,z; active Hamilton quaternion. Reject nonunit input rather
// than interpreting arbitrary four-vectors as orientations.
Rotation rotation_from_quaternion(double w,double x,double y,double z);

struct PhaseGridView
{
    const double* u_edges=nullptr;
    Vec3 ki{};
    Frame frame{};
    std::uint32_t nt=0,np=0;
};
inline constexpr double pole_roundoff_tolerance=256*DBL_EPSILON;
struct PhaseBin
{
    std::uint64_t index=0;
    double mu=0;
    int pole=0; // +1 forward, -1 backward, 0 ordinary
    bool valid=false,boundary_snapped=false;
};
ICE_HD inline PhaseBin phase_bin(PhaseGridView g,Vec3 direction) noexcept
{
    PhaseBin b{};
    if(!g.u_edges||!g.nt||!g.np||!is_unit(direction))return b;
    const double x=g.frame.e0.dot(direction),y=g.frame.e1.dot(direction),z=g.ki.dot(direction);
    const double n=::sqrt(x*x+y*y+z*z);
    if(!(n>0)||!finite_value(n)||::fabs(n-1)>512*DBL_EPSILON)return b;
    const double mu=z/n,s=::hypot(x,y)/n;
    b.mu=mu;
    if(s<=pole_roundoff_tolerance)
    {
        b.pole=mu>=0?1:-1;b.valid=true;return b;
    }
    // Stable at the two poles. The SAVED binary64 u edges define membership.
    double u=mu>=0?(0.5*s*s)/(1+mu):1-(0.5*s*s)/(1-mu);
    if(!(u>=0&&u<=1))return b;
    std::uint32_t lo=0,hi=g.nt;
    while(lo+1<hi)
    {
        const auto mid=lo+(hi-lo)/2;
        if(u>=g.u_edges[mid])lo=mid;else hi=mid;
    }
    // Canonical half-open boundary assignment within roundoff only. The cap
    // prevents this from becoming a finite-width filter on extremely fine grids.
    for(unsigned candidate=0;candidate<2;++candidate)
    {
        const auto e=lo+candidate;
        if(e==0||e>=g.nt)continue;
        const double width=::fmin(g.u_edges[e]-g.u_edges[e-1],g.u_edges[e+1]-g.u_edges[e]);
        const double tolerance=::fmin(64*DBL_EPSILON*::fmax(::fabs(u),::fabs(g.u_edges[e])),1e-7*width);
        if(::fabs(u-g.u_edges[e])<=tolerance)
        {b.boundary_snapped=(u!=g.u_edges[e]);lo=e;break;}
    }
    double phi=::atan2(y,x);
    if(phi>=pi)phi-=2*pi;
    const double v=(phi+pi)/(2*pi);
    double index=v*g.np;
    const double nearest=::floor(index+0.5);
    if(::fabs(index-nearest)<=::fmin(64*DBL_EPSILON*g.np,1e-7))
    {b.boundary_snapped=b.boundary_snapped||(index!=nearest);index=nearest;}
    auto j=static_cast<std::uint32_t>(index);
    if(j==g.np)j=0; // periodic seam, including the last rounding ulp
    if(j>=g.np)return b;
    b.index=std::uint64_t(lo)*g.np+j;b.valid=true;return b;
}
struct HistogramStats
{
    double point_mass=0,point_axial=0,forward_mass=0,backward_mass=0;
    double absorbed_addend_mass=0; // diagnostic, NOT an error bound
    std::uint64_t samples=0,invalid_samples=0,absorbed_addends=0,boundary_snapped_samples=0;
};
struct HistogramMoments {double mass=0,axial=0;};
struct PhaseAccumulateParams
{
    PhaseGridView grid{};
    Rotation rotation{};
    const OutgoingSample* outgoing=nullptr;
    double* cell_mass=nullptr;
    HistogramStats* statistics=nullptr;
    std::uint64_t count=0;
    double area_weight_per_ray=0;
};
struct PhaseFinishParams
{
    PhaseGridView grid{};
    const double* cell_mass=nullptr;
    const HistogramStats* statistics=nullptr;
    double* density=nullptr;
    HistogramMoments* partials=nullptr;
};
static_assert(std::is_trivially_copyable_v<PhaseAccumulateParams>);
static_assert(std::is_trivially_copyable_v<PhaseFinishParams>);
} // namespace iceCrystal
