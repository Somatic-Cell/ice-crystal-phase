#pragma once

#include <rainbow/host_device.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

namespace rainbow
{


struct Complex32
{
    float real = 0.0f;
    float imag = 0.0f;

    [[nodiscard]]
    HOST_DEVICE constexpr Complex32 operator+(const Complex32 rhs) const noexcept
    {
        return {real + rhs.real, imag + rhs.imag};
    }
    
    [[nodiscard]]
    HOST_DEVICE constexpr Complex32 operator-(const Complex32 rhs) const noexcept
    {
        return {real - rhs.real, imag - rhs.imag};
    }

    [[nodiscard]]
    HOST_DEVICE Complex32 operator*(const Complex32 rhs) const noexcept
    {
        // (a + ib)(c + id) = (ac - bd) + i(ad + bc)
        // 各成分と一方の積と加算を fmaf で一回の丸め処理で済ませる
        return {
            ::fmaf(real, rhs.real, - (imag * rhs.imag)),    // ac - bd
            ::fmaf(real, rhs.imag, imag * rhs.real)         // ad + bc
        };
    }

    [[nodiscard]]
    HOST_DEVICE constexpr Complex32 operator*(const float scale) const noexcept
    {
        return {real * scale, imag * scale};
    }

    HOST_DEVICE constexpr Complex32& operator+=(const Complex32 rhs) noexcept
    {
        real += rhs.real;
        imag += rhs.imag;
        return *this;
    }

    HOST_DEVICE constexpr Complex32& operator-=(const Complex32 rhs) noexcept
    {
        real -= rhs.real;
        imag -= rhs.imag;
        return *this;
    }
    
    HOST_DEVICE Complex32& operator*=(const Complex32 rhs) noexcept
    {
        const Complex32 product = (*this) * rhs;
        real = product.real;
        imag = product.imag;
        return *this;
    }

    [[nodiscard]]
    HOST_DEVICE constexpr Complex32 conjugate() const noexcept
    {
        return {real, -imag};
    }

    [[nodiscard]]
    HOST_DEVICE float squared_norm() const noexcept
    {
        return fmaf(real, real, imag * imag);
    }


};

static_assert(sizeof(float) == 4);
static_assert(std::numeric_limits<float>::is_iec559);
static_assert(std::numeric_limits<float>::digits == 24);
static_assert(std::is_standard_layout_v<Complex32>);
static_assert(std::is_trivially_copyable_v<Complex32>);
static_assert(sizeof(Complex32) == 8);
static_assert(alignof(Complex32) == 4);
static_assert(offsetof(Complex32, real) == 0);
static_assert(offsetof(Complex32, imag) == 4);

}