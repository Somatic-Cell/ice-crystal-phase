#pragma once

#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <rainbow/raindrop_shape.hpp>
#include <rainbow/phase_cycles.hpp>
#include <rainbow/field32.hpp>

namespace rainbow
{

enum TraceDiagnostic : std::uint32_t { NoDiagnostic=0, IntersectionFallback=1, SnellFallback=2 };

enum class PathFamily : std::uint32_t { R=0, TT=1, TRT=2, TRRT=3 };
enum class VertexStatus : std::uint32_t
{
    Unwritten=0, Valid=1, Miss=2, TotalInternalReflection=3,
    UnresolvedIntersection=4, InvalidInterface=5, PhaseOverflow=6
};
struct OutgoingVertex
{
    Vec3 position_drop{};       // 無次元．world/mm へ移さず後段へ渡す．
    Vec3 direction_drop{};
    Vec3 basis_x{};             // basis_y = direction_drop x basis_x．
    Field32 field{};            // 伝播・焦線位相，patch の面積比は未適用．
    PhaseCycles optical_cycles{};
    VertexStatus status=VertexStatus::Unwritten;
    std::uint32_t diagnostics=0;
};
static_assert(sizeof(OutgoingVertex)==68 && alignof(OutgoingVertex)==4);
static_assert(std::is_trivially_copyable_v<OutgoingVertex>);
static_assert(std::is_standard_layout_v<OutgoingVertex>);

// shape / 格子配置は全 vertex 共通．一つの入射 Jones 状態を決定論的に追跡する．
struct RaindropTraceConfig
{
    RaindropShape shape{};
    Vec3 incident_direction{1,0,0};
    Vec3 incident_basis_x{0,0,1};
    Field32 incident_field{{1,0},{0,0}};
    float radius_mm=0.4f;
    float wavelength_nm=700.0f;
    float exterior_index=1.0f;
    float interior_index=1.3314f;
    float grid_half_extent=1.01f;
    float reference_distance=2.0f;
    float outgoing_reference_distance=2.0f;
    std::uint32_t grid_width=129;
    std::uint32_t grid_height=129;
    HOST_DEVICE std::uint32_t vertex_count() const noexcept {return grid_width*grid_height;}
};
static_assert(sizeof(RaindropTraceConfig)==116);
static_assert(std::is_trivially_copyable_v<RaindropTraceConfig>);
static_assert(std::is_standard_layout_v<RaindropTraceConfig>);

// OptixTraversableHandle と同じ幅の値を保持するだけ．数値ヘッダは OptiX に依存させない．
struct alignas(8) RaindropTraceParams
{
    // constant 変数自体は trivial default construction にする．
    // default member initializers を持つ config は別の GPU buffer に置く．
    OutgoingVertex* vertices;
    std::uint64_t traversable;
    const RaindropTraceConfig* config;
};
static_assert(sizeof(void*)==8, "This project requires a 64-bit build.");
static_assert(sizeof(RaindropTraceParams)==24);
static_assert(alignof(RaindropTraceParams)==8);
static_assert(std::is_trivially_default_constructible_v<RaindropTraceParams>);
static_assert(std::is_standard_layout_v<RaindropTraceParams>);
static_assert(offsetof(RaindropTraceParams,config)==16);
static_assert(std::is_trivially_copyable_v<RaindropTraceParams>);
}
