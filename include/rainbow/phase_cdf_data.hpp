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
// Additional storage filtering is explicitly separate from the optical
// diffraction approximation. Zero dimensions inherit the query grid.
struct PhaseStorageSettings
{
    std::uint32_t theta_count = 0, phi_count = 0;
    double gaussian_sigma_degrees = 0;
    double gaussian_support_sigma = 4;
    double maximum_coarsening_tv = 1; // structural default, NOT accuracy certification
};
struct PhaseStorageStatistics
{
    double g_source = 0, g_filtered = 0, g_stored = 0;
    double gaussian_integral_relative_change = 0;
    double aggregation_integral_relative_change = 0;
    double coarsening_tv = 0, cdf_mass = 0;
    bool gaussian_underresolved = false;
};
struct PhaseMomentSum
{
    double mass, axial, tv;
    std::uint64_t weight_underflow;
};
struct PhaseStorageParams
{
    PhaseDensityView input{};
    const double* fine_edges = nullptr;
    const double* fine_values = nullptr;
    const double* coarse_values = nullptr;
    const double* phi_cdf = nullptr;
    const double* theta_cdf = nullptr;
    const double* coarse_edges = nullptr;
    double* output = nullptr;
    PhaseMomentSum* moment_partials = nullptr;
    std::uint32_t nt = 0, np = 0, out_nt = 0, out_np = 0;
    double divisor = 1, sigma_rad = 0, support_sigma = 4;
    double fine_mass = 1, coarse_mass = 1;
};
static_assert(std::is_trivially_copyable_v<PhaseStorageParams>);
static_assert(std::is_trivially_default_constructible_v<PhaseMomentSum>);
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
