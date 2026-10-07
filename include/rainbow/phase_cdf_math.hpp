#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
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

struct PhaseCdfSample
{
    double u = 0, v = 0, pdf_omega = 0;
    std::uint32_t theta_cell = 0, phi_cell = 0;
};

// Caller validates complete tables once. Each call checks only its random input
// and selected CDF interval; no PDF reconstruction or renormalization.
[[nodiscard]] inline PhaseCdfSample sample_phase_cdf(
    std::span<const double> phi, std::span<const double> theta,
    std::span<const double> edges, double xi_phi, double xi_theta)
{
    if(phi.size() < 2 || edges.size() < 2 ||
       theta.size() != (phi.size()-1)*edges.size() ||
       !(xi_phi >= 0 && xi_phi < 1 && xi_theta >= 0 && xi_theta < 1))
        throw std::invalid_argument("Invalid CDF layout/random coordinates.");
    const auto np = phi.size()-1;
    const auto nt = edges.size()-1;
    const auto j = std::size_t(std::upper_bound(phi.begin(), phi.end(), xi_phi)-phi.begin()-1);
    if(j >= np) throw std::runtime_error("Invalid marginal CDF.");
    const auto column = theta.subspan(j*(nt+1), nt+1);
    const auto i = std::size_t(std::upper_bound(column.begin(), column.end(), xi_theta)-column.begin()-1);
    if(i >= nt) throw std::runtime_error("Invalid conditional CDF.");
    const double a = phi[j+1]-phi[j], b = column[i+1]-column[i];
    const double du = edges[i+1]-edges[i];
    if(!(a > 0 && b > 0 && du > 0)) throw std::runtime_error("Selected an empty interval.");
    const double r0 = (xi_phi-phi[j])/a, r1 = (xi_theta-column[i])/b;
    double u = edges[i] + r1*du;
    double v = (double(j)+r0)/double(np);
    // Keep rounded values in the selected half-open cell. Not CDF clamping.
    u = (std::min)(u, std::nextafter(edges[i+1], edges[i]));
    v = (std::min)(v, std::nextafter(double(j+1)/double(np), double(j)/double(np)));
    return {u, v, (a*b*double(np))/(4*phase_cdf_pi*du),
            static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(j)};
}
} // namespace rainbow
