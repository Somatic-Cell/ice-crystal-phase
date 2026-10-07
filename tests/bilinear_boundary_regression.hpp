#pragma once

#include <rainbow/bilinear_patch_intersection.hpp>
#include <bit>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace rainbow::tests
{
namespace boundary_regression_detail
{
inline void require(bool condition, const char* message)
{
    if(!condition) throw std::runtime_error(message);
}
inline Vec3 from_bits(std::uint32_t x, std::uint32_t y, std::uint32_t z)
{
    return {std::bit_cast<float>(x), std::bit_cast<float>(y), std::bit_cast<float>(z)};
}
inline BilinearPatchGeometry rectangle(float x0, float x1, float y0, float y1)
{
    return {{{x0, y0, 1.0f}, {x1, y0, 1.0f},
             {x0, y1, 1.0f}, {x1, y1, 1.0f}}};
}
inline void require_one_regular_hit(const BilinearRayHits& result, const char* message)
{
    require(result.count == 1u, message);
    require((result.flags & (BilinearUnresolved | BilinearInvalidInput)) == 0u, message);
    require((result.hits[0].flags & (BilinearBoundary | BilinearSingular)) == 0u, message);
}

// CPU reconstruction of the R patch at grid=3001, cell=(2132,1330),
// radius=1 mm, inclination=20 degrees, query=(623,796) on a 900x1800 grid.
// This is NOT a capture of the user's GPU buffers. Store the input bits to
// isolate the query regression from unrelated future changes in ray tracing.
inline void test_guard_equality()
{
    const BilinearPatchGeometry patch{{
        from_bits(0xbedf3e39u, 0x3ef02275u, 0x3f449f68u),
        from_bits(0xbedeb0fau, 0x3eefdebeu, 0x3f44dc13u),
        from_bits(0xbedfac65u, 0x3eefa235u, 0x3f44a739u),
        from_bits(0xbedf1f2du, 0x3eef5e8fu, 0x3f44e3e4u)}};
    const Vec3 direction = from_bits(0xbedf00f7u, 0x3eef81b8u, 0x3f44e1c1u);
    const auto result = BilinearPatchIntersector::intersect(patch, direction);
    require_one_regular_hit(result, "Interior R root was misclassified at the FP32 guard equality.");
    // Independent long-double triple-product reference (GCC/x86-64).
    // The tolerances below cover storage of u/v in binary32, not a new
    // acceptance tolerance in the production intersector.
    require(std::abs(double(result.hits[0].u) - 0.999995492741538133) < 1e-7,
            "R witness u differs from the reference root.");
    require(std::abs(double(result.hits[0].v) - 0.725773698683431621) < 1e-7,
            "R witness v differs from the reference root.");
}

inline void test_packed_endpoint()
{
    constexpr float extent = 0x1p-10f;
    constexpr float tiny = 0x1p-37f;
    const auto patch = rectangle(-extent, tiny, -extent, extent);
    const auto result = BilinearPatchIntersector::intersect(patch, {0, 0, 1});
    require_one_regular_hit(result, "Packing an interior u to 1.0f created a false boundary.");
    // Exact root: u = 2^27/(2^27+1), v = 1/2. It is inside the square,
    // but correctly rounded binary32 u equals 1. Do not perturb the ray.
    const double reference_u = 134217728.0 / 134217729.0;
    require(reference_u < 1.0 && static_cast<float>(reference_u) == 1.0f,
            "Synthetic endpoint-rounding witness is not representable as expected.");
    require(result.hits[0].u == 1.0f && result.hits[0].v == 0.5f,
            "Packed endpoint witness changed.");
    require((result.hits[0].flags & BilinearUsedFp64) != 0u,
            "Near-edge witness did not use the existing FP64 fallback.");

    const auto neighbor = rectangle(tiny, extent, -extent, extent);
    const auto adjacent = BilinearPatchIntersector::intersect(neighbor, {0, 0, 1});
    require(adjacent.count == 0u && (adjacent.flags & BilinearUnresolved) == 0u,
            "The adjacent patch must not acquire the interior root.");

    const auto patch_v = rectangle(-extent, extent, -extent, tiny);
    const auto result_v = BilinearPatchIntersector::intersect(patch_v, {0, 0, 1});
    require_one_regular_hit(result_v, "Packing an interior v to 1.0f created a false boundary.");
    require(result_v.hits[0].u == 0.5f && result_v.hits[0].v == 1.0f,
            "Packed v endpoint witness changed.");
}

inline void test_real_boundaries_and_misses()
{
    constexpr float extent = 0x1p-10f;
    constexpr float tiny = 0x1p-37f;
    const BilinearPatchGeometry edges[] = {
        rectangle(-extent, 0, -extent, extent),
        rectangle(0, extent, -extent, extent),
        rectangle(-extent, extent, -extent, 0),
        rectangle(-extent, extent, 0, extent)};
    for(const auto& patch : edges)
    {
        const auto result = BilinearPatchIntersector::intersect(patch, {0, 0, 1});
        require(result.count == 1u, "True closed edge was lost.");
        require((result.flags & BilinearUnresolved) == 0u, "True planar edge became unresolved.");
        require((result.hits[0].flags & BilinearBoundary) != 0u,
                "True edge must remain marked for ownership handling.");
    }
    const float endpoints[2][2] = {{-extent, 0}, {0, extent}};
    for(const auto& x : endpoints) for(const auto& y : endpoints)
    {
        const auto patch = rectangle(x[0], x[1], y[0], y[1]);
        const auto result = BilinearPatchIntersector::intersect(patch, {0, 0, 1});
        require(result.count == 1u && (result.hits[0].flags & BilinearBoundary) != 0u,
                "True corner must remain marked for ownership handling.");
    }
    const BilinearPatchGeometry outside[] = {
        rectangle(-extent, -tiny, -extent, extent),
        rectangle(tiny, extent, -extent, extent),
        rectangle(-extent, extent, -extent, -tiny),
        rectangle(-extent, extent, tiny, extent)};
    for(const auto& patch : outside)
    {
        const auto result = BilinearPatchIntersector::intersect(patch, {0, 0, 1});
        require(result.count == 0u && (result.flags & BilinearUnresolved) == 0u,
                "Outside root was incorrectly accepted.");
    }
}
} // namespace boundary_regression_detail

inline void test_bilinear_boundary_regressions()
{
    boundary_regression_detail::test_guard_equality();
    boundary_regression_detail::test_packed_endpoint();
    boundary_regression_detail::test_real_boundaries_and_misses();
}
} // namespace rainbow::tests
