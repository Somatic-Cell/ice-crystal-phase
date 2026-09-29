#pragma once

#include <rainbow/patch_query_data.hpp>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace rainbow
{
// constant-memory object 自体は trivial default construction にする．
// offsets[i] から summaries[i].hit_count 件が有効（容量との差は padding）．
struct alignas(8) PatchQueryLaunchParams
{
    std::uint64_t traversable;
    const OutgoingVertex* vertices;
    const OutgoingPatch* patches;
    const Vec3* directions;
    const std::uint64_t* offsets;
    PatchQuerySummary* summaries;
    PatchQueryHit* hits;
    std::uint32_t patch_count;
    std::uint32_t vertex_count;
    std::uint32_t direction_count;
    std::uint32_t write_pass;
};
static_assert(sizeof(PatchQueryLaunchParams)==72 && alignof(PatchQueryLaunchParams)==8);
static_assert(offsetof(PatchQueryLaunchParams,patch_count)==56);
static_assert(std::is_trivially_default_constructible_v<PatchQueryLaunchParams>);
static_assert(std::is_standard_layout_v<PatchQueryLaunchParams>);
static_assert(std::is_trivially_copyable_v<PatchQueryLaunchParams>);
} // namespace rainbow
