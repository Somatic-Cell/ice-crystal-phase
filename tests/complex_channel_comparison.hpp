#pragma once

#include <algorithm>
#include <cmath>

namespace unpolarized_tests
{
// CPU/GPU comparison of one complex output-polarization channel.
// Compare the complex DIFFERENCE, not the difference of magnitudes: a phase
// error must remain observable. Check s and p separately at the call site.
// This helper changes only the test metric, never any simulator value.
[[nodiscard]]
inline bool near_complex_channel(
    const double actual_real, const double actual_imag,
    const double reference_real, const double reference_imag,
    const double relative, const double absolute = 1e-12) noexcept
{
    if(!std::isfinite(actual_real) || !std::isfinite(actual_imag)
       || !std::isfinite(reference_real) || !std::isfinite(reference_imag)
       || !std::isfinite(relative) || !std::isfinite(absolute)
       || relative < 0.0 || absolute < 0.0)
    {
        return false;
    }
    if(actual_real == reference_real && actual_imag == reference_imag)
    {
        return true;
    }

    const double error = std::hypot(
        actual_real - reference_real, actual_imag - reference_imag);
    const double scale = (std::max)(
        std::hypot(actual_real, actual_imag),
        std::hypot(reference_real, reference_imag));
    const double allowed = absolute + relative * scale;
    // An overflow in the comparison itself is not evidence of agreement.
    return std::isfinite(error) && std::isfinite(scale)
        && std::isfinite(allowed) && error <= allowed;
}
} // namespace unpolarized_tests
