#pragma once
#include <rainbow/patch_optics_data.hpp>
#include <cstdint>
#include <type_traits>

namespace rainbow
{
// Companion buffers. ABI v2 includes the second incident response column.
struct FocalPhaseConfig
{
    std::uint32_t grid_width = 0, grid_height = 0;
    float grid_half_extent = 0;
    std::uint32_t reserved = 0;
    // Paper: positive derivative adds one crossing relative to negative.
    // Absolute offsets between families are NOT specified there. All-zero is
    // an explicit project convention, NOT an exact caustic-crossing count.
    std::uint32_t quarter_turn_offsets[4] = {0, 0, 0, 0};
};
enum FocalOpticalFlags : std::uint32_t
{
    FocalNone = 0, FocalInputError = 1u, FocalInputPending = 2u,
    FocalDerivativePending = 4u, FocalInvalidGeometry = 8u,
    FocalArithmeticError = 16u
};
inline constexpr std::uint32_t focal_error_mask = FocalInputError | FocalInvalidGeometry | FocalArithmeticError;
inline constexpr std::uint32_t focal_pending_mask = FocalInputPending | FocalDerivativePending;
struct FocalOpticalResult
{
    OpticalField64 field{};
    double intensity_s = 0, intensity_p = 0;
    double family_incoherent[4] = {};
    std::uint32_t family_hits[4] = {};
    std::uint32_t corrected_hits = 0, extra_quarter_turn_hits = 0;
    std::uint32_t flags = 0, first_problem_patch_id = 0xffffffffu;
    // field is input column 0; field_second is column 1. Intensities are the
    // incoherent input-state average in unpolarized mode, not |field|^2 alone.
    OpticalField64 field_second{};
    [[nodiscard]] HOST_DEVICE bool valid() const noexcept { return flags == 0u; }
};
struct WaveOpticsParams
{
    PatchOpticsParams optical{};
    FocalPhaseConfig focal{};
    FocalOpticalResult* focal_results = nullptr;
};
struct DiffractionConfig
{
    double primary_sigma_rad = 0;
    double support_sigma = 4;
    double transition_contrast = 0.25;
    double inner_blend_sigma = 3;
};
enum RainbowTransitionFlags : std::uint32_t
{
    TransitionNone = 0, TransitionPrimary = 1u, TransitionSecondary = 2u,
    TransitionUnknown = 4u
};
// Edge between this theta row and the next, at fixed phi; no edge on last row.
struct RainbowTransition { std::uint32_t flags = 0; };
enum DiffractionFlags : std::uint32_t
{
    DiffractionNone = 0, DiffractionFiltered = 1u, DiffractionUnderresolved = 2u,
    DiffractionInputUnavailable = 4u, DiffractionStencilUnavailable = 8u,
    DiffractionDetectionUnavailable = 16u, DiffractionInvalidData = 32u
};
inline constexpr std::uint32_t diffraction_unavailable_mask = DiffractionInputUnavailable | DiffractionStencilUnavailable
    | DiffractionDetectionUnavailable | DiffractionInvalidData;
struct DiffractionResult
{
    double intensity_s = 0, intensity_p = 0, sigma_rad = 0, blend = 0;
    std::uint32_t flags = 0, transition_kind = 0;
    [[nodiscard]] HOST_DEVICE bool valid() const noexcept
    { return (flags & diffraction_unavailable_mask) == 0u; }
};
struct DiffractionParams
{
    const FocalOpticalResult* focal = nullptr;
    RainbowTransition* transitions = nullptr;
    DiffractionResult* results = nullptr;
    std::uint32_t theta_count = 0, phi_count = 0;
    DiffractionConfig config{};
};
struct WaveOpticsSettings
{
    std::uint32_t focal_quarter_turn_offsets[4] = {0, 0, 0, 0};
    // Zero selects Table II. Positive = explicit user-specified extension.
    double primary_sigma_degrees = 0;
    double transition_contrast = 0.25;
};
struct WaveOpticsStatistics
{
    std::uint32_t focal_errors = 0, focal_pending = 0;
    std::uint32_t primary_transitions = 0, secondary_transitions = 0, unknown_transitions = 0;
    std::uint32_t filtered_directions = 0, unavailable_directions = 0, underresolved_directions = 0;
    std::uint64_t corrected_hits = 0, extra_quarter_turn_hits = 0;
};
static_assert(sizeof(FocalPhaseConfig) == 32);
static_assert(sizeof(FocalOpticalResult) == 144);
static_assert(sizeof(WaveOpticsParams) == 160);
static_assert(sizeof(DiffractionResult) == 40);
static_assert(sizeof(DiffractionParams) == 64);
static_assert(std::is_trivially_copyable_v<WaveOpticsParams>);
static_assert(std::is_trivially_copyable_v<FocalOpticalResult>);
static_assert(std::is_trivially_copyable_v<DiffractionParams>);
static_assert(std::is_standard_layout_v<WaveOpticsParams>);
} // namespace rainbow
