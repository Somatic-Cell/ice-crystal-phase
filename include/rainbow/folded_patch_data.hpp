#pragma once

#include <rainbow/wave_optics_data.hpp>
#include <cstdint>
#include <type_traits>

namespace rainbow
{
// Nonregular geometry is NOT relabelled as regular. These companion records
// describe a separate, explicit finite-branch-density approximation.
enum FoldedPatchFlags : std::uint32_t
{
    FoldedReady = 0,
    FoldedInvalidInput = 1u << 0,
    FoldedHemisphereUncertified = 1u << 1,
    FoldedNoSignChange = 1u << 2,
    FoldedTopologyUnresolved = 1u << 3,
    FoldedQuadratureUnconverged = 1u << 4,
    FoldedInvalidMetric = 1u << 5
};
struct FoldedPatchConfig
{
    // Estimated quadrature accuracy, NOT a rigorous error bound on optics.
    double relative_tolerance = 1e-8;
    std::uint32_t maximum_panels = 2048;
    std::uint32_t maximum_depth = 20;
};
struct FoldedBranch
{
    // The branch is the indicated connected component of sign*J>0.
    // Each vertical section is a single interval. Its u projection is [u0,u1].
    double u0 = 0, u1 = 0;
    double area_fraction = 0, solid_angle_sr = 0;
    double representative_u = 0, representative_v = 0;
    std::int32_t sign = 0;
    std::uint32_t focal_flags = FocalInputPending;
    std::uint32_t quarter_turns = 0, extra_quarter_turn = 0;
};
struct FoldedPatchRecord
{
    std::uint32_t compact_index = 0, patch_id = 0;
    std::uint32_t flags = FoldedInvalidInput, branch_count = 0;
    // Four values, order 00,10,01,11. Scaled by max(abs(J_corner)).
    double jacobian[4] = {};
    FoldedBranch branches[4] = {};
    double largest_estimated_relative_error = 0;
    std::uint32_t panels = 0, reserved = 0;
};
struct FoldedPatchView
{
    const std::uint32_t* record_indices = nullptr; // by compact patch index
    const FoldedPatchRecord* records = nullptr;
    std::uint32_t patch_count = 0, record_count = 0;
    [[nodiscard]] HOST_DEVICE bool enabled() const noexcept
    { return record_indices != nullptr && records != nullptr; }
    [[nodiscard]] HOST_DEVICE const FoldedPatchRecord* find(std::uint32_t compact) const noexcept
    {
        if(!enabled() || compact >= patch_count) return nullptr;
        const auto i = record_indices[compact];
        if(i >= record_count) return nullptr;
        const auto* r = records + i;
        return r->compact_index == compact ? r : nullptr;
    }
};
// The wave parameters include the optional second incident-response column.
struct FoldedPrepareParams
{
    const OutgoingVertex* vertices = nullptr;
    const OutgoingPatch* patches = nullptr;
    std::uint32_t* record_indices = nullptr;
    FoldedPatchRecord* records = nullptr;
    std::uint32_t* written_count = nullptr;
    std::uint32_t vertex_count = 0, patch_count = 0, capacity = 0;
    Vec3 incident_direction{};
    std::uint32_t with_focal = 0;
    FocalPhaseConfig focal{};
    FoldedPatchConfig config{};
};
struct FoldedOpticsParams
{
    WaveOpticsParams wave{};
    FoldedPatchView folded{};
    std::uint32_t with_focal = 0, reserved = 0;
};
struct FoldedPatchStatistics
{
    std::uint32_t candidates = 0, prepared = 0, unresolved = 0;
    std::uint32_t branch_count = 0, focal_pending_branches = 0;
    std::uint64_t evaluated_hits = 0;
    double largest_estimated_relative_error = 0;
};
static_assert(std::is_trivially_copyable_v<FoldedPatchRecord>);
static_assert(std::is_standard_layout_v<FoldedPatchRecord>);
static_assert(std::is_trivially_copyable_v<FoldedPrepareParams>);
static_assert(std::is_trivially_copyable_v<FoldedOpticsParams>);
static_assert(sizeof(FoldedPatchView) == 24);
static_assert(sizeof(FoldedBranch) == 64 && sizeof(FoldedPatchRecord) == 320);
static_assert(sizeof(FoldedPrepareParams) == 120 && sizeof(FoldedOpticsParams) == 192);
} // namespace rainbow
