#pragma once

#include <rainbow/host_device.hpp>
#include <cmath>
#include <type_traits>

namespace rainbow
{

// 正規化された粒子局所座標の小さな値型．double は交差の曖昧な場合と参照計算用．
template<class T> struct Vec3T
{
    static_assert(std::is_same_v<T,float> || std::is_same_v<T,double>);
    T x = T(0), y = T(0), z = T(0);

    HOST_DEVICE Vec3T operator+(const Vec3T b) const noexcept { return {x+b.x,y+b.y,z+b.z}; }
    HOST_DEVICE Vec3T operator-(const Vec3T b) const noexcept { return {x-b.x,y-b.y,z-b.z}; }
    HOST_DEVICE Vec3T operator-() const noexcept { return {-x,-y,-z}; }
    HOST_DEVICE Vec3T operator*(const T s) const noexcept { return {x*s,y*s,z*s}; }
    HOST_DEVICE Vec3T operator/(const T s) const noexcept { return {x/s,y/s,z/s}; }
    HOST_DEVICE T dot(const Vec3T b) const noexcept
    {
        if constexpr(std::is_same_v<T,float>) return ::fmaf(x,b.x,::fmaf(y,b.y,z*b.z));
        else return ::fma(x,b.x,::fma(y,b.y,z*b.z));
    }
    HOST_DEVICE Vec3T cross(const Vec3T b) const noexcept
    {
        // FMA 一方だけでは v.cross(v) が厳密な零にならないため，ここは対称な差．
        return {y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x};
    }
    HOST_DEVICE T length() const noexcept
    {
        if constexpr(std::is_same_v<T,float>) return ::sqrtf(dot(*this));
        else return ::sqrt(dot(*this));
    }
    HOST_DEVICE Vec3T normalized() const noexcept { return *this / length(); }
    HOST_DEVICE bool is_finite() const noexcept
    {
        // host-only std::numeric_limits の関数を device 側から呼ばない．
        T limit;
        if constexpr(std::is_same_v<T,float>) limit=T(0x1.fffffep127f);
        else limit=T(0x1.fffffffffffffp1023);
        return x >= -limit && x <= limit && y >= -limit && y <= limit && z >= -limit && z <= limit;
    }
    HOST_DEVICE Vec3T at(const Vec3T direction, const T t) const noexcept
    {
        if constexpr(std::is_same_v<T,float>) return {::fmaf(t,direction.x,x),::fmaf(t,direction.y,y),::fmaf(t,direction.z,z)};
        else return {::fma(t,direction.x,x),::fma(t,direction.y,y),::fma(t,direction.z,z)};
    }
    template<class U> HOST_DEVICE Vec3T<U> cast() const noexcept
    { return {static_cast<U>(x),static_cast<U>(y),static_cast<U>(z)}; }
};
using Vec3 = Vec3T<float>;

// e0, e1, w が右手系になる transverse frame．normal incidence の基底選択にも使う．
struct TransverseFrame
{
    Vec3 e0;
    Vec3 e1;
    HOST_DEVICE static TransverseFrame from_direction(const Vec3 w) noexcept
    {
        const Vec3 axis = ::fabsf(w.z) < 0.9f ? Vec3{0,0,1} : Vec3{1,0,0};
        const Vec3 e0 = (axis-w*w.dot(axis)).normalized();
        return {e0,w.cross(e0).normalized()};
    }
};
static_assert(sizeof(Vec3)==12 && alignof(Vec3)==4);
static_assert(std::is_trivially_copyable_v<Vec3> && std::is_standard_layout_v<Vec3>);
}
