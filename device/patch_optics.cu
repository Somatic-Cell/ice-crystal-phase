#include <rainbow/patch_optical_evaluator.hpp>
#include <cstdint>

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
