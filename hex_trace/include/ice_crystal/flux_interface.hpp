#pragma once
#include <ice_crystal/math.hpp>

namespace iceCrystal
{
struct FluxRay { Vec3 direction; Frame frame; JonesFlux response; };
struct FluxBranches
{
    FluxRay reflected,transmitted;
    double balance_error=0;
    bool valid=false,has_transmission=false,total_internal_reflection=false;
};
ICE_HD inline Complex tir_phase(double a,double b) noexcept
{
    const double s=max_value(::fabs(a),::fabs(b));
    a/=s;b/=s;
    const double den=a*a+b*b;
    return {(a*a-b*b)/den,-2*a*b/den};
}

// Nonabsorbing, nonmagnetic, ISOTROPIC dielectrics. The normal points into the
// incident medium. A JonesFlux norm represents POWER, not electric intensity:
// t_flux = sqrt(n_t*cos_t/(n_i*cos_i))*t_E.
// The equivalent symmetric formula below avoids a divergent grazing prefactor.
// Complex TIR phases are retained even though different paths are NOT interfered.
ICE_HD inline FluxBranches split_flux(const FluxRay& ray,Vec3 n,double ni,double nt) noexcept
{
    FluxBranches b{};
    if(!ray.frame.valid(ray.direction)||!is_unit(n)||!ray.response.is_finite()||
        !(ni>0&&nt>0&&finite_value(ni)&&finite_value(nt)))return b;
    double ci=-ray.direction.dot(n);
    if(ci<0 || ci>1+128*DBL_EPSILON)return b;
    if(ci>1)ci=1; // unit-vector roundoff only, never an optical coefficient repair
    if(ni==nt)
    {
        b.transmitted=ray;b.has_transmission=true;
        b.reflected.direction=ray.direction+n*(2*ci);
        if(!normalize(b.reflected.direction,b.reflected.direction))return b;
        if(!make_frame(b.reflected.direction,b.reflected.frame))return b;
        b.valid=true;return b;
    }
    Vec3 s=ray.direction.cross(n);
    const double scale=s.max_abs();
    const Vec3 scaled=scale>0?s/scale:Vec3{};
    const double sine=scale*::sqrt(scaled.dot(scaled));
    if(sine>0) { if(!normalize(s,s))return b; }
    else s=ray.frame.e0;
    Vec3 pin{};
    if(!normalize(ray.direction.cross(s),pin))return b;
    const auto j=change_frame(ray.response,ray.frame,{s,pin});
    const double ns=max_value(ni,nt),a=ni/ns,c=nt/ns,eta=a/c;
    if(!finite_value(eta)||!(a>0&&c>0))return b;
    const double st=eta*sine;
    if(!finite_value(st))return b;
    if(!normalize(ray.direction+n*(2*ci),b.reflected.direction))return b;
    b.reflected.frame.e0=s;
    if(!normalize(b.reflected.direction.cross(s),b.reflected.frame.e1))return b;
    Complex rs{},rp{};
    if(st>1)
    {
        const double inverse=1/st;
        const double beta=st*::sqrt((1-inverse)*(1+inverse));
        rs=tir_phase(a*ci,c*beta);rp=tir_phase(c*ci,a*beta);
        b.total_internal_reflection=true;
    }
    else
    {
        const double ct=::sqrt((1-st)*(1+st));
        const double ds=a*ci+c*ct,dp=c*ci+a*ct;
        if(!(ds>0&&dp>0))return b;
        rs={(a*ci-c*ct)/ds,0};rp={(c*ci-a*ct)/dp,0};
        const double numerator=2*::sqrt(a*ci)*::sqrt(c*ct);
        const double ts=numerator/ds,tp=numerator/dp;
        b.transmitted.response=scale_rows(j,{ts,0},{tp,0});
        if(!normalize((ray.direction+n*ci)*eta-n*ct,b.transmitted.direction))return b;
        b.transmitted.frame.e0=s;
        if(!normalize(b.transmitted.direction.cross(s),b.transmitted.frame.e1))return b;
        b.has_transmission=true;
    }
    b.reflected.response=scale_rows(j,rs,rp);
    const double before=ray.response.power();
    const double after=b.reflected.response.power()+(b.has_transmission?b.transmitted.response.power():0);
    b.balance_error=after-before;
    b.valid=b.reflected.frame.valid(b.reflected.direction)&&b.reflected.response.is_finite()&&
        (!b.has_transmission || (b.transmitted.frame.valid(b.transmitted.direction)&&b.transmitted.response.is_finite()))&&
        finite_value(b.balance_error) && ::fabs(b.balance_error)<=2048*DBL_EPSILON*max_value(before,DBL_MIN);
    return b;
}
} // namespace iceCrystal
