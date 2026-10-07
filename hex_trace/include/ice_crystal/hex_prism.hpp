#pragma once
#include <ice_crystal/math.hpp>

namespace iceCrystal
{
// Right regular hexagonal prism, centered at zero, longitudinal axis +y.
// Side face j (0..5): n=(cos(j*pi/3),0,sin(j*pi/3)), n.x <= sqrt(3)*R/2.
// Face 6: +y; face 7: -y. R is circumradius, L is FULL length, both mm.
// Geometry calculations use x/R so uniform rescaling does not change tolerances.
struct HexPrism
{
    double circumradius_mm=0, length_mm=0;
    ICE_HD double half_length() const noexcept { return 0.5*(length_mm/circumradius_mm); }
    ICE_HD bool valid() const noexcept
    {
        return finite_value(circumradius_mm)&&finite_value(length_mm)&&circumradius_mm>0&&length_mm>0&&
            finite_value(half_length())&&half_length()>0 &&
            finite_value(circumradius_mm*circumradius_mm)&&circumradius_mm*circumradius_mm>0 &&
            finite_value(circumradius_mm*length_mm)&&circumradius_mm*length_mm>0;
    }
    ICE_HD Vec3 normal(unsigned f) const noexcept
    {
        switch(f)
        {
        case 0:return {1,0,0}; case 1:return {0.5,0,sqrt3/2};
        case 2:return {-0.5,0,sqrt3/2};case 3:return {-1,0,0};
        case 4:return {-0.5,0,-sqrt3/2};case 5:return {0.5,0,-sqrt3/2};
        case 6:return {0,1,0};case 7:return {0,-1,0};default:return {};
        }
    }
    ICE_HD Vec3 vertex(unsigned f) const noexcept
    {
        switch(f%6)
        {
        case 0:return {sqrt3/2,0,0.5};case 1:return {0,0,1};
        case 2:return {-sqrt3/2,0,0.5};case 3:return {-sqrt3/2,0,-0.5};
        case 4:return {0,0,-1}; default:return {sqrt3/2,0,-0.5};
        }
    }
    ICE_HD double offset(unsigned f) const noexcept { return f<6?sqrt3/2:half_length(); }
    ICE_HD double face_area_mm2(unsigned f) const noexcept
    { return f<6?circumradius_mm*length_mm:(3*sqrt3/2)*circumradius_mm*circumradius_mm; }
    ICE_HD double surface_area_mm2() const noexcept
    { return 6*face_area_mm2(0)+2*face_area_mm2(6); }
    ICE_HD double projected_area_mm2(Vec3 k) const noexcept
    {
        Sum a{};
        for(unsigned f=0;f<8;++f)
            a.add(face_area_mm2(f)*max_value(0,-normal(f).dot(k)));
        return a.value;
    }
};

ICE_HD inline double geometry_tolerance(const HexPrism& h,Vec3 p,double t=0) noexcept
{
    return 128*DBL_EPSILON*max_value(1,max_value(h.half_length(),max_value(p.max_abs(),::fabs(t))));
}
struct EntrySample
{
    Vec3 position; // normalized body coordinates x/R; NOT mm
    unsigned face=8;
    bool valid=false;
};
ICE_HD inline EntrySample sample_entry(const HexPrism& h,Vec3 k,double u0,double u1,double u2) noexcept
{
    EntrySample s{};
    if(!h.valid()||!is_unit(k)||!(u0>=0&&u0<1&&u1>=0&&u1<1&&u2>=0&&u2<1)) return s;
    const double a=h.projected_area_mm2(k);
    if(!(a>0)||!finite_value(a)) return s;
    const double target=u0*a;
    double accumulated=0;
    unsigned last=8;
    for(unsigned f=0;f<8;++f)
    {
        const double w=h.face_area_mm2(f)*max_value(0,-h.normal(f).dot(k));
        if(!(w>0)) continue;
        last=f; accumulated+=w;
        if(target<accumulated) { s.face=f; break; }
    }
    // Only the rounding-sized final interval; never select a back-facing face.
    if(s.face==8) s.face=last;
    if(s.face==8) return s;
    if(s.face<6)
    {
        const auto n=h.normal(s.face);
        const Vec3 tangent{-n.z,0,n.x};
        s.position=n*(sqrt3/2)+tangent*(u1-0.5)+Vec3{0,(2*u2-1)*h.half_length(),0};
    }
    else
    {
        const double sector=6*u1;
        const auto j=static_cast<unsigned>(sector);
        const double r=::sqrt(sector-double(j));
        s.position=(h.vertex(j)*(1-u2)+h.vertex(j+1)*u2)*r;
        s.position.y=(s.face==6?1:-1)*h.half_length();
    }
    s.valid=s.position.is_finite();return s;
}

enum class HitStatus : std::uint32_t { Regular, Miss, AmbiguousBoundary, Invalid };
struct BoundaryHit
{
    Vec3 position;
    double distance=0; // distance / R; direction must be unit length
    unsigned face=8;
    HitStatus status=HitStatus::Invalid;
};
ICE_HD inline bool strict_face_point(const HexPrism& h,Vec3 p,unsigned face) noexcept
{
    if(face>=8 || !p.is_finite()) return false;
    const double tol=geometry_tolerance(h,p);
    for(unsigned f=0;f<8;++f)
    {
        const double gap=h.offset(f)-h.normal(f).dot(p);
        if(f==face) { if(::fabs(gap)>tol) return false; }
        else if(!(gap>tol)) return false;
    }
    return true;
}

// For an internal ray starting on a KNOWN face. No origin shift or epsilon
// t-min. An unresolved edge/vertex is reported, not assigned an arbitrary normal.
ICE_HD inline BoundaryHit next_boundary(const HexPrism& h,Vec3 p,Vec3 d,unsigned start_face) noexcept
{
    BoundaryHit hit{};
    if(!h.valid()||!is_unit(d)||start_face>=8) return hit;
    if(!strict_face_point(h,p,start_face)) { hit.status=HitStatus::AmbiguousBoundary;return hit; }
    if(h.normal(start_face).dot(d)>128*DBL_EPSILON) return hit;
    double best=DBL_MAX,second=DBL_MAX;
    for(unsigned f=0;f<8;++f)
    {
        if(f==start_face) continue;
        const double cosine=h.normal(f).dot(d);
        if(!(cosine>0)) continue;
        const double t=(h.offset(f)-h.normal(f).dot(p))/cosine;
        if(!(t>0)||!finite_value(t)) continue;
        if(t<best) {second=best;best=t;hit.face=f;}
        else if(t<second)second=t;
    }
    if(hit.face==8)return hit;
    hit.distance=best;hit.position=p+d*best;
    const double tol=geometry_tolerance(h,hit.position,best);
    if(second-best<=tol || !strict_face_point(h,hit.position,hit.face))
        hit.status=HitStatus::AmbiguousBoundary;
    else hit.status=HitStatus::Regular;
    return hit;
}

// External half-space clipping. Returns the first nonnegative boundary. Input
// and output are normalized body coordinates; useful for independent ray tests.
ICE_HD inline BoundaryHit intersect_prism(const HexPrism& h,Vec3 p,Vec3 d) noexcept
{
    BoundaryHit hit{};
    if(!h.valid()||!p.is_finite()||!is_unit(d))return hit;
    double near_t=-DBL_MAX,far_t=DBL_MAX;
    unsigned near_f=8,far_f=8;
    for(unsigned f=0;f<8;++f)
    {
        const double gap=h.offset(f)-h.normal(f).dot(p),den=h.normal(f).dot(d);
        if(den==0)
        {
            if(gap<0){hit.status=HitStatus::Miss;return hit;}
            continue;
        }
        const double t=gap/den;
        if(den<0 && t>near_t){near_t=t;near_f=f;}
        if(den>0 && t<far_t){far_t=t;far_f=f;}
    }
    if(far_t<near_t || far_t<0){hit.status=HitStatus::Miss;return hit;}
    hit.distance=near_t>=0?near_t:far_t;hit.face=near_t>=0?near_f:far_f;
    if(hit.face==8||!finite_value(hit.distance))return hit;
    hit.position=p+d*hit.distance;
    hit.status=strict_face_point(h,hit.position,hit.face)?HitStatus::Regular:HitStatus::AmbiguousBoundary;
    return hit;
}
} // namespace iceCrystal
