#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace rainbow
{
inline constexpr double phase_cdf_pi = 3.141592653589793238462643383279502884;

// Store these exact binary64 coordinates, and use their adjacent differences
// in BOTH construction and sampling. This avoids cross-language cos rounding.
[[nodiscard]] inline std::vector<double> make_phase_u_edges(std::uint32_t n)
{
    if(n == 0) throw std::invalid_argument("theta_count must be positive.");
    std::vector<double> edges(std::size_t(n) + 1u);
    edges.front() = 0;
    edges.back() = 1;
    for(std::uint32_t i = 1; i < n; ++i)
    {
        const auto k = (std::min)(i, n - i);
        const double s = std::sin((0.5 * phase_cdf_pi) * double(k) / double(n));
        edges[i] = i <= n / 2u ? s*s : 1.0 - s*s;
    }
    for(std::uint32_t i = 0; i < n; ++i)
        if(!(edges[i+1] > edges[i]))
            throw std::range_error("Angular cells cannot be represented in binary64.");
    return edges;
}

} // namespace rainbow
