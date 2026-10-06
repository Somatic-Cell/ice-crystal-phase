#pragma once

// Host-only diagnostics. Never used by the production classification/solver.
#include <rainbow/patch_geometry.hpp>
#include <array>
#include <cmath>
#include <limits>

namespace rainbow
{
// Re-evaluate the SAME stored FP32 corners. FP64 does not recover precision
// lost during tracing; the exact offline audit uses these corners as inputs.
struct PatchRegularityAudit
{
    struct Bound { double lo = 0, hi = 0; };
    struct Result
    {
        std::array<Bound,4> jacobian32{}, jacobian64{}, hemisphere32{}, hemisphere64{};
        std::array<double,4> jacobian_value64{};
        Vec3 hemisphere_axis{}; // exactly the FP32 candidate used by the original code
        int host_orientation32 = 0, interval_orientation64 = 0;
        float host_signed_omega32 = 0;
        double signed_omega64 = 0;
        bool finite_input = false;
    };

    [[nodiscard]] static Result inspect(const BilinearPatchGeometry& geometry) noexcept
    {
        Result r{};
        for(const auto p : geometry.corners) if(!p.is_finite()) return r;
        r.finite_input = true;
        r.hemisphere_axis=(geometry.corners[0]+geometry.corners[1])
                          +(geometry.corners[2]+geometry.corners[3]);
        // Replay of the original FP32 predicate, on the host, NOT a GPU trace.
        r.host_orientation32=geometry.certified_orientation();
        r.host_signed_omega32=geometry.signed_solid_angle();
        evaluate_intervals<float>(geometry,r.hemisphere_axis,r.jacobian32,r.hemisphere32);
        evaluate_intervals<double>(geometry,r.hemisphere_axis,r.jacobian64,r.hemisphere64);
        bool positive=true,negative=true,hemisphere=true;
        for(unsigned k=0;k<4;++k)
        {
            positive &= r.jacobian64[k].lo>0;
            negative &= r.jacobian64[k].hi<0;
            hemisphere &= r.hemisphere64[k].lo>0;
        }
        r.interval_orientation64=hemisphere ? (positive?1:(negative?-1:0)) : 0;
        const auto a=geometry.corners[0].cast<double>(), b=geometry.corners[1].cast<double>();
        const auto c=geometry.corners[2].cast<double>(), d=geometry.corners[3].cast<double>();
        r.jacobian_value64={a.dot((b-a).cross(c-a)),b.dot((b-a).cross(d-b)),
                            c.dot((d-c).cross(c-a)),d.dot((d-c).cross(d-b))};
        r.signed_omega64=triangle_omega(a,b,d)+triangle_omega(a,d,c);
        return r;
    }

    [[nodiscard]] static double jacobian_at(const BilinearPatchGeometry& g,
                                            const double u,const double v) noexcept
    {
        const auto a=g.corners[0].cast<double>(), b=g.corners[1].cast<double>();
        const auto c=g.corners[2].cast<double>(), d=g.corners[3].cast<double>();
        const auto p=(a*(1-u)+b*u)*(1-v)+(c*(1-u)+d*u)*v;
        const auto pu=(b-a)*(1-v)+(d-c)*v;
        const auto pv=(c-a)*(1-u)+(d-b)*u;
        return p.dot(pu.cross(pv));
    }

private:
    // One representable step outward after EVERY elementary operation.
    // Must compile without fast-math, reassociation, implicit FMA or FTZ/DAZ.
    template<class T> struct Interval
    {
        T lo,hi;
        static Interval point(const T x) noexcept {return {x,x};}
        static T down(const T x) noexcept {return std::nextafter(x,-std::numeric_limits<T>::infinity());}
        static T up(const T x) noexcept {return std::nextafter(x,std::numeric_limits<T>::infinity());}
        Interval operator+(const Interval q) const noexcept {return {down(lo+q.lo),up(hi+q.hi)};}
        Interval operator-(const Interval q) const noexcept {return {down(lo-q.hi),up(hi-q.lo)};}
        Interval operator*(const Interval q) const noexcept
        {
            const T x=lo*q.lo,y=lo*q.hi,z=hi*q.lo,w=hi*q.hi;
            return {down(std::fmin(std::fmin(x,y),std::fmin(z,w))),
                    up(std::fmax(std::fmax(x,y),std::fmax(z,w)))};
        }
    };
    template<class T> struct Vector
    {
        Interval<T> x,y,z;
        static Vector point(const Vec3 p) noexcept
        {return {Interval<T>::point(T(p.x)),Interval<T>::point(T(p.y)),Interval<T>::point(T(p.z))};}
        Vector operator-(const Vector q) const noexcept {return {x-q.x,y-q.y,z-q.z};}
        Vector cross(const Vector q) const noexcept
        {return {y*q.z-z*q.y,z*q.x-x*q.z,x*q.y-y*q.x};}
        Interval<T> dot(const Vector q) const noexcept {return x*q.x+y*q.y+z*q.z;}
    };
    template<class T> static void evaluate_intervals(const BilinearPatchGeometry& g,const Vec3 axis,
        std::array<Bound,4>& j,std::array<Bound,4>& h) noexcept
    {
        using V=Vector<T>;
        const V a=V::point(g.corners[0]),b=V::point(g.corners[1]);
        const V c=V::point(g.corners[2]),d=V::point(g.corners[3]);
        const Interval<T> q[4]={a.dot((b-a).cross(c-a)),b.dot((b-a).cross(d-b)),
                                c.dot((d-c).cross(c-a)),d.dot((d-c).cross(d-b))};
        const V n=V::point(axis);
        const Interval<T> dots[4]={n.dot(a),n.dot(b),n.dot(c),n.dot(d)};
        for(unsigned k=0;k<4;++k)
        {j[k]={double(q[k].lo),double(q[k].hi)}; h[k]={double(dots[k].lo),double(dots[k].hi)};}
    }
    static double triangle_omega(const Vec3T<double> a,const Vec3T<double> b,
                                 const Vec3T<double> c) noexcept
    {
        const double la=a.length(),lb=b.length(),lc=c.length();
        const double denominator=la*lb*lc+a.dot(b)*lc+b.dot(c)*la+c.dot(a)*lb;
        return 2*std::atan2(a.dot((b-a).cross(c-a)),denominator);
    }
};
} // namespace rainbow
