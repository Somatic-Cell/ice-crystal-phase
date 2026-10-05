#include <rainbow/patch_optical_evaluator.hpp>
#include <cstdint>
#include <rainbow/rainbow_diffraction.hpp>

// Ordinary CUDA, not an OptiX program. One direction has one owning thread.
// No fixed hit limit; the evaluator consumes summary.hit_count, not capacity.
extern "C" __global__
void evaluate_patch_optics(const rainbow::PatchOpticsParams params)
{
    const std::uint64_t index = std::uint64_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if(index < params.direction_count)
        params.results[index] = rainbow::PatchOpticalEvaluator::evaluate_direction(
            params, static_cast<std::uint32_t>(index));
}

// Same per-hit evaluator, one traversal of the existing hit list, two outputs.
extern "C" __global__
void evaluate_focal_optics(const rainbow::WaveOpticsParams params)
{
    const std::uint64_t i = std::uint64_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if(i < params.optical.direction_count)
        params.optical.results[i] = rainbow::PatchOpticalEvaluator::evaluate_direction(
            params.optical, static_cast<std::uint32_t>(i), &params.focal, params.focal_results + i);
}

extern "C" __global__
void detect_rainbow_transitions(const rainbow::DiffractionParams params)
{
    const std::uint64_t i = std::uint64_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if(i < std::uint64_t(params.theta_count) * params.phi_count)
        params.transitions[i] = rainbow::RainbowDiffraction::detect(params, static_cast<std::uint32_t>(i));
}

extern "C" __global__
void filter_rainbow_diffraction(const rainbow::DiffractionParams params)
{
    const std::uint64_t i = std::uint64_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if(i < std::uint64_t(params.theta_count) * params.phi_count)
        params.results[i] = rainbow::RainbowDiffraction::filter(params, static_cast<std::uint32_t>(i));
}
