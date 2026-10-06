#pragma once

#include <rainbow/field32.hpp>
#include <cstdint>
#include <type_traits>

namespace rainbow
{
// A response operator, NOT a Jones vector for unpolarized light.
// Column j is the output field for the j-th unit incident Jones basis state.
// All vertices use the SAME input basis. Rows use each ray's transverse frame.
struct JonesResponse32
{
    Field32 column[2]{};

    [[nodiscard]] HOST_DEVICE static constexpr JonesResponse32 identity() noexcept
    {
        JonesResponse32 result{};
        result.column[0].x = {1.0f, 0.0f};
        result.column[1].y = {1.0f, 0.0f};
        return result;
    }

    [[nodiscard]] HOST_DEVICE Field32 apply_to(const Field32 input) const noexcept
    {
        return {column[0].x * input.x + column[1].x * input.y,
                column[0].y * input.x + column[1].y * input.y};
    }

    [[nodiscard]] HOST_DEVICE bool is_finite() const noexcept
    {
        for(unsigned j = 0; j < 2; ++j)
        {
            const auto& f = column[j];
            constexpr float limit = 0x1.fffffep127f;
            if(!(f.x.real >= -limit && f.x.real <= limit
                 && f.x.imag >= -limit && f.x.imag <= limit
                 && f.y.real >= -limit && f.y.real <= limit
                 && f.y.imag >= -limit && f.y.imag <= limit)) return false;
        }
        return true;
    }
};

// Explicit representation tag. An absent second column MUST NOT silently turn
// an unpolarized calculation into a single-input calculation.
enum class IncidentPolarization : std::uint32_t
{
    SingleJones = 0, // retained for reference tests / explicit library use
    Unpolarized = 1
};

static_assert(sizeof(JonesResponse32) == 32 && alignof(JonesResponse32) == 4);
static_assert(std::is_standard_layout_v<JonesResponse32>);
static_assert(std::is_trivially_copyable_v<JonesResponse32>);
} // namespace rainbow
