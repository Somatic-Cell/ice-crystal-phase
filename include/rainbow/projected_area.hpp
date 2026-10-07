#pragma once

#include <array>
#include <cstdint>

namespace rainbow
{
// CPU-only geometry; uses the SAME dimensionless radial coefficients as the
// optical solver, promoted to binary64. Polar axis is particle-space -y.
struct ProjectedAreaOptions
{
    double relative_tolerance = 1e-10;
    double absolute_tolerance_dimensionless = 1e-12;
    std::uint32_t maximum_quadrature_panels = 32768;
};

struct ProjectedAreaResult
{
    double area_mm2 = 0;
    double estimated_absolute_error_mm2 = 0;
    double estimated_relative_error = 0;
    std::uint32_t quadrature_evaluations = 0;
    std::uint32_t convexity_intervals = 0;
    bool analytic_sphere = false;
};

// r(chi)/a = 1 + sum_{n=0}^7 c[n] cos(n chi).
// The surface-projection identity requires convexity. It is checked with
// outward-rounded interval arithmetic over the WHOLE polar domain, not by
// point sampling. Nonconvex/unverified shapes throw; no convex-hull substitute.
// The quadrature error is an ESTIMATE, not a rigorous enclosure of the area.
// No wavelength, ray grid, angular grid, or CDF participates in this calculation.
[[nodiscard]] ProjectedAreaResult compute_projected_area(
    const std::array<double, 8>& coefficients,
    double radius_mm,
    const std::array<double, 3>& incident_propagation,
    const ProjectedAreaOptions& options = {});
}
