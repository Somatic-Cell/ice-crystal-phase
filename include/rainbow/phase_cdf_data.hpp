#pragma once

#include <rainbow/patch_optics_data.hpp>
#include <rainbow/wave_optics_data.hpp>
#include <cstdint>
#include <type_traits>

namespace rainbow
{
enum class PhaseDensityStage : std::uint32_t
{
    Scalar = 0, Incoherent = 1, Path = 2, Focal = 3, Diffraction = 4
};

// Optional scalar input supports independent numerical tests and other solvers.
// Product code uses optical/focal/diffraction with their original validity flags.
struct PhaseDensityView
{
    const double* scalar = nullptr;
    const PatchOpticalResult* optical = nullptr;
    const FocalOpticalResult* focal = nullptr;
    const DiffractionResult* diffraction = nullptr;
    std::uint64_t count = 0;
    PhaseDensityStage stage = PhaseDensityStage::Scalar;
    std::uint32_t reserved = 0;
};

// Trivial default construction is required for CUDA shared-memory arrays.
struct PhaseCdfReduction
{
    double maximum;
    double l1_error, lost_mass, maximum_cell_error;
    std::uint64_t invalid_values, incomplete, underresolved;
    std::uint64_t malformed, lost_cells, weight_underflow;
    std::uint64_t first_problem;
};
struct PhaseCdfPolicy
{
    double maximum_l1_error = 1e-10;
    double maximum_lost_mass = 1e-12;
    bool allow_underresolved = false;
};
struct PhaseCdfBuildParams
{
    PhaseDensityView input{};
    const double* u_edges = nullptr;
    const PhaseCdfReduction* input_summary = nullptr;
    double* phi_cdf = nullptr;
    double* theta_cdf = nullptr;
    double* column_sums = nullptr;
    double* total = nullptr;
    PhaseCdfReduction* partials = nullptr;
    PhaseCdfReduction* output_summary = nullptr;
    std::uint32_t theta_count = 0, phi_count = 0;
};
struct TraceValidationParams
{
    const OutgoingVertex* vertices = nullptr;
    const Field32* second = nullptr;
    std::uint64_t count = 0;
    PhaseCdfReduction* partials = nullptr;
};
static_assert(sizeof(PhaseDensityView) == 48);
static_assert(sizeof(PhaseCdfReduction) == 88);
static_assert(std::is_trivially_default_constructible_v<PhaseCdfReduction>);
static_assert(std::is_trivially_copyable_v<PhaseCdfBuildParams>);
static_assert(std::is_standard_layout_v<PhaseCdfBuildParams>);
} // namespace rainbow
