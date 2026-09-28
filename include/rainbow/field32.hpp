#pragma once

#include <rainbow/complex32.hpp>

#include <cstddef>
#include <type_traits>

namespace rainbow
{

struct Field32
{
    Complex32 x = {};
    Complex32 y = {};

    [[nodiscard]]
    HOST_DEVICE constexpr Field32 operator+(const Field32 rhs) const noexcept
    {
        return {x + rhs.x, y + rhs.y};
    }
    
    [[nodiscard]]
    HOST_DEVICE constexpr Field32 operator-(const Field32 rhs) const noexcept
    {
        return {x - rhs.x, y - rhs.y};
    }

    HOST_DEVICE constexpr Field32 operator+=(const Field32 rhs) noexcept
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }
    
    HOST_DEVICE constexpr Field32 operator-=(const Field32 rhs) noexcept
    {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }
    
    [[nodiscard]]
    HOST_DEVICE Field32 operator*(const float scale) const noexcept
    {
        return {x *scale, y * scale};
    }
    
    HOST_DEVICE Field32& operator*=(const Complex32 factor) noexcept
    {
        x *= factor;
        y *= factor;
        return *this;
    }

    [[nodiscard]]
    HOST_DEVICE float squared_norm() const noexcept
    {
        return x.squared_norm() + y.squared_norm();
    }

    [[nodiscard]]
    HOST_DEVICE Field32 rotated_basis(
        const float cosine,
        const float sine
    ) const noexcept
    {
        return {
            {
                ::fmaf(cosine, x.real, sine * y.real),
                ::fmaf(cosine, x.imag, sine * y.imag)
            },
            {
                ::fmaf(-sine, x.real, cosine * y.real),
                ::fmaf(-sine, x.imag, cosine * y.imag)
            }
        };
    }
};

static_assert(std::is_standard_layout_v<Field32>);
static_assert(std::is_trivially_copyable_v<Field32>);
static_assert(sizeof(Field32) == 16);
static_assert(alignof(Field32) == 4);
static_assert(offsetof(Field32, x) == 0);
static_assert(offsetof(Field32, y) == 8);

} // namespace rainbow