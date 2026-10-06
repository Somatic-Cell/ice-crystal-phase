#pragma once

#include <rainbow/patch_query_data.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace rainbow
{

// This stage deliberately omits focal-line phase and diffraction.  These flags
// describe evaluation of the supplied hits, NOT completeness of the phase LUT.
enum PatchOpticalFlags : std::uint32_t
{
    PatchOpticalNone              = 0,
    PatchOpticalQueryError        = 1u << 0,
    PatchOpticalInvalidInput      = 1u << 1,
    PatchOpticalFrameFailure      = 1u << 2,
    PatchOpticalPhaseFailure      = 1u << 3,
    PatchOpticalArithmeticFailure = 1u << 4,
    PatchOpticalRefinementPending = 1u << 5,
    PatchOpticalBoundaryPending   = 1u << 6,
    PatchOpticalSingularPending   = 1u << 7
};
inline constexpr std::uint32_t patch_optical_error_mask =
    PatchOpticalQueryError | PatchOpticalInvalidInput | PatchOpticalFrameFailure
    | PatchOpticalPhaseFailure | PatchOpticalArithmeticFailure;
inline constexpr std::uint32_t patch_optical_pending_mask =
    PatchOpticalRefinementPending | PatchOpticalBoundaryPending | PatchOpticalSingularPending;

// Namespace-scope scalar constant: no host-only standard-library call in a kernel.
inline constexpr double patch_optical_nan = std::numeric_limits<double>::quiet_NaN();

// s = perpendicular to the scattering plane; p = outgoing_direction cross s.
// This is a value/accumulation type, not a replacement for the stored Field32.
struct OpticalField64
{
    double s_real = 0.0, s_imag = 0.0;
    double p_real = 0.0, p_imag = 0.0;
};

struct PatchOpticalResult
{
    // Historical field names retained for source/ABI compatibility. With the
    // optional folded-patch view, these sums also include resolved branch hits.
    // A numeric error makes all values NaN; unsupported hits remain explicit.
    OpticalField64 regular_partial_path_field{};
    double regular_partial_incoherent_s = 0.0;
    double regular_partial_incoherent_p = 0.0;
    double regular_partial_path_s = 0.0;
    double regular_partial_path_p = 0.0;

    std::uint32_t hit_count = 0;
    std::uint32_t evaluated_hits = 0;
    std::uint32_t rejected_hits = 0;
    std::uint32_t refinement_hits = 0;
    std::uint32_t boundary_hits = 0;
    std::uint32_t singular_hits = 0;
    std::uint32_t flags = PatchOpticalNone;
    std::uint32_t query_flags = 0;
    std::uint32_t first_problem_patch_id = 0xffffffffu;
    // Formerly reserved ABI word. Tests that use it as an output guard still work.
    std::uint32_t reserved = 0;
    [[nodiscard]] HOST_DEVICE std::uint32_t folded_evaluated_hits() const noexcept
    { return reserved; }

    [[nodiscard]] HOST_DEVICE bool known_hits_complete() const noexcept
    {
        // Does NOT certify missing source cells, angular/grid convergence,
        // focal-line phase, diffraction, or agreement with Mie theory.
        return (flags & (patch_optical_error_mask | patch_optical_pending_mask)) == 0u
            && rejected_hits == 0u && evaluated_hits == hit_count;
    }
};

// Passed BY VALUE to the ordinary CUDA kernel. All arrays must belong to the
// same CUDA context and remain unchanged until that kernel finishes.
// Offset differences are capacities; only summaries[i].hit_count slots are read.
struct PatchOpticsParams
{
    const OutgoingVertex* vertices = nullptr;
    const OutgoingPatch* patches = nullptr;
    const Vec3* directions = nullptr;
    const std::uint64_t* offsets = nullptr; // direction_count + 1 entries
    const PatchQuerySummary* summaries = nullptr;
    const PatchQueryHit* hits = nullptr;
    PatchOpticalResult* results = nullptr;
    std::uint64_t hit_storage_count = 0;
    std::uint32_t vertex_count = 0;
    std::uint32_t patch_count = 0;
    std::uint32_t direction_count = 0;
    std::uint32_t reserved = 0;
    Vec3 incident_direction{};
    Vec3 incident_basis_x{};
};

static_assert(sizeof(OpticalField64) == 32);
static_assert(sizeof(PatchOpticalResult) == 104 && alignof(PatchOpticalResult) == 8);
static_assert(offsetof(PatchOpticalResult, hit_count) == 64);
static_assert(sizeof(void*) == 8);
static_assert(sizeof(PatchOpticsParams) == 104 && alignof(PatchOpticsParams) == 8);
static_assert(offsetof(PatchOpticsParams, incident_direction) == 80);
static_assert(std::is_standard_layout_v<PatchOpticalResult>);
static_assert(std::is_trivially_copyable_v<PatchOpticalResult>);
static_assert(std::is_standard_layout_v<PatchOpticsParams>);
static_assert(std::is_trivially_copyable_v<PatchOpticsParams>);

} // namespace rainbow
