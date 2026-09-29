#pragma once

#include <rainbow/raindrop_trace_data.hpp>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace rainbow
{

// One status per logical input cell and path family. Missing cells are retained
// in this array; compaction must never erase the evidence of a boundary.
enum class PatchCellStatus : std::uint32_t
{
    Unwritten = 0,
    RegularPositive = 1,
    RegularNegative = 2,
    NeedsRefinement = 3,       // Four rays exist, but a regular spherical map was not certified.
    MissingCorners = 4,       // One to three outgoing rays exist.
    NoOutgoingCorners = 5,    // Corner evidence only; NOT a proof that the whole cell is empty.
    InvalidVertex = 6,        // Unwritten / unresolved / invalid interface / phase overflow.
    InvalidGeometry = 7      // A supposedly valid vertex contains unusable data.
};

inline constexpr unsigned int patch_status_count = 8;
inline constexpr unsigned int patch_threads_per_block = 256;

[[nodiscard]]
HOST_DEVICE constexpr bool has_patch_geometry(const PatchCellStatus status) noexcept
{
    return status == PatchCellStatus::RegularPositive
        || status == PatchCellStatus::RegularNegative
        || status == PatchCellStatus::NeedsRefinement;
}

// CUDA has no dependency on OptiX headers. src/patch_accel.cpp checks every field
// against OptixAabb. The allocation (not necessarily the C++ type) is 8-byte aligned.
struct PatchAabb
{
    float min_x, min_y, min_z;
    float max_x, max_y, max_z;
};

struct OutgoingPatch
{
    // Tensor-product order: (00,10,01,11), NOT boundary order (00,10,11,01).
    // These are global indices into the borrowed, path-major OutgoingVertex array.
    std::uint32_t vertex_indices[4] = {};
    std::uint32_t patch_id = 0; // family*cells_per_path + iy*(grid_width-1) + ix.
    PatchCellStatus status = PatchCellStatus::Unwritten;
    float signed_solid_angle_sr = 0.0f;
    float incident_area_drop2 = 0.0f;

    [[nodiscard]]
    HOST_DEVICE constexpr bool has_regular_spherical_map() const noexcept
    {
        return status == PatchCellStatus::RegularPositive
            || status == PatchCellStatus::RegularNegative;
    }
};

// Kernel parameter, passed BY VALUE. No CUDA/OptiX handles, no owning pointers.
struct PatchBuildLayout
{
    std::uint32_t grid_width = 0;
    std::uint32_t grid_height = 0;
    std::uint32_t vertices_per_path = 0;
    std::uint32_t cells_per_path = 0;
    std::uint32_t cell_count = 0; // All four path families.
    float incident_area_drop2 = 0.0f;

    // This stage uses exactly the existing square emitting grid.
    // It does not change the silhouette bounds or divide energy by active cells.
    [[nodiscard]]
    static bool try_make(const RaindropTraceConfig& config, PatchBuildLayout& result) noexcept
    {
        if(config.grid_width < 2 || config.grid_height < 2
           || !(config.grid_half_extent > 0.0f)
           || !(config.grid_half_extent <= (std::numeric_limits<float>::max)()))
            return false;
        const std::uint64_t vertices = std::uint64_t(config.grid_width)*config.grid_height;
        const std::uint64_t cells = std::uint64_t(config.grid_width-1)*(config.grid_height-1);
        if(vertices > (std::numeric_limits<std::uint32_t>::max)()/4u
           || cells > (std::numeric_limits<std::uint32_t>::max)()/4u)
            return false;
        const float du = 2.0f*config.grid_half_extent/float(config.grid_width-1);
        const float dv = 2.0f*config.grid_half_extent/float(config.grid_height-1);
        const float area = du*dv;
        if(!(area > 0.0f && area <= (std::numeric_limits<float>::max)())) return false;
        result = {config.grid_width, config.grid_height,
                  static_cast<std::uint32_t>(vertices), static_cast<std::uint32_t>(cells),
                  static_cast<std::uint32_t>(4*cells), area};
        return true;
    }

    HOST_DEVICE void corner_indices(const std::uint32_t id, std::uint32_t (&out)[4]) const noexcept
    {
        const std::uint32_t family = id/cells_per_path;
        const std::uint32_t cell = id%cells_per_path;
        const std::uint32_t ix = cell%(grid_width-1);
        const std::uint32_t iy = cell/(grid_width-1);
        const std::uint32_t first = family*vertices_per_path + iy*grid_width + ix;
        out[0]=first; out[1]=first+1; out[2]=first+grid_width; out[3]=first+grid_width+1;
    }
};

struct PatchBlockSummary
{
    std::uint32_t counts[patch_status_count] = {};
};

static_assert(sizeof(PatchAabb) == 24 && alignof(PatchAabb) == 4);
static_assert(sizeof(OutgoingPatch) == 32 && alignof(OutgoingPatch) == 4);
static_assert(sizeof(PatchBuildLayout) == 24);
static_assert(sizeof(PatchBlockSummary) == 32);
static_assert(std::is_trivially_copyable_v<OutgoingPatch> && std::is_standard_layout_v<OutgoingPatch>);
static_assert(std::is_trivially_copyable_v<PatchAabb> && std::is_standard_layout_v<PatchAabb>);
static_assert(std::is_trivially_copyable_v<PatchBuildLayout> && std::is_standard_layout_v<PatchBuildLayout>);
static_assert(std::is_trivially_copyable_v<PatchBlockSummary>);
static_assert(offsetof(OutgoingPatch, patch_id) == 16);
static_assert(offsetof(OutgoingPatch, signed_solid_angle_sr) == 24);
static_assert(offsetof(PatchAabb, max_x) == 12);

} // namespace rainbow
