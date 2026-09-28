#pragma once

#include <cstdint>
#include <cmath>
#include <rainbow/phase_cycles.hpp>
#include <rainbow/vec3.hpp>

namespace rainbow
{

// FP32 の hi+lo で中間積・除算の余りを保持する．絶対 OPL を double で保存しない．
// 正規化粒子座標，正の可視波長，mm サイズの入力範囲を対象とする．
struct FloatPair
{
    float hi=0.0f,lo=0.0f;
    HOST_DEVICE static FloatPair sum(const float a,const float b) noexcept
    {
        const float s=a+b,v=s-a;
        return {s,(a-(s-v))+(b-v)};
    }
    HOST_DEVICE static FloatPair product(const float a,const float b) noexcept
    {
        const float p=a*b;
        return {p,::fmaf(a,b,-p)};
    }
    HOST_DEVICE FloatPair operator+(const FloatPair b) const noexcept
    {
        const auto s=sum(hi,b.hi);
        const auto t=sum(lo,b.lo);
        const auto v=sum(s.lo,t.hi);
        const auto w=sum(s.hi,v.hi);
        return sum(w.hi,((w.lo+v.lo)+t.lo));
    }
    HOST_DEVICE FloatPair operator-() const noexcept {return {-hi,-lo};}
    HOST_DEVICE FloatPair operator-(const FloatPair b) const noexcept {return *this+(-b);}
    HOST_DEVICE FloatPair operator*(const FloatPair b) const noexcept
    {
        const auto p=product(hi,b.hi);
        return sum(p.hi,::fmaf(hi,b.lo,::fmaf(lo,b.hi,p.lo))+lo*b.lo);
    }
    HOST_DEVICE FloatPair divided_by(const float denominator) const noexcept
    {
        const float q=hi/denominator;
        const float residual=::fmaf(-q,denominator,hi)+lo;
        return sum(q,residual/denominator);
    }
    HOST_DEVICE FloatPair square_root() const noexcept
    {
        const float r=::sqrtf(hi);
        if(r==0.0f) return {};
        return sum(r,(::fmaf(-r,r,hi)+lo)/(2.0f*r));
    }
    HOST_DEVICE static FloatPair dot(const Vec3 a,const Vec3 b) noexcept
    {
        return product(a.x,b.x)+product(a.y,b.y)+product(a.z,b.z);
    }
    HOST_DEVICE static FloatPair distance(const Vec3 a,const Vec3 b) noexcept
    {
        const auto x=sum(a.x,-b.x),y=sum(a.y,-b.y),z=sum(a.z,-b.z);
        return (x*x+y*y+z*z).square_root();
    }
};

// q = (radius_mm * 10^6 / wavelength_nm) * n * distance_drop．
// 単一 float の巨大な q から fraction を取り出すことはしない．
struct OpticalPathAccumulator
{
    PhaseCycles phase{};
    FloatPair cycles_per_drop_unit{};

    HOST_DEVICE static OpticalPathAccumulator make(const float radius_mm,const float wavelength_nm) noexcept
    {
        OpticalPathAccumulator p{};
        p.cycles_per_drop_unit=FloatPair::product(radius_mm,1000000.0f).divided_by(wavelength_nm);
        return p;
    }
    HOST_DEVICE bool add_distance(const FloatPair distance,const float refractive_index) noexcept
    {
        const FloatPair q=(cycles_per_drop_unit*FloatPair{refractive_index,0})*distance;
        if(!(q.hi>=0.0f && q.hi<0x1p31f && ::fabsf(q.lo)<0x1p31f)) return false;
        const float integral=::floorf(q.hi);
        PhaseCycles increment{},correction{},candidate=phase;
        if(!increment.try_set_parts(static_cast<std::int32_t>(integral),q.hi-integral)
           ||!correction.try_set_parts(0,q.lo)
           ||!increment.try_add(correction)||!candidate.try_add(increment)) return false;
        phase=candidate;
        return true;
    }
    HOST_DEVICE bool add_segment(const Vec3 from,const Vec3 to,const float n) noexcept
    { return add_distance(FloatPair::distance(from,to),n); }
};
}
