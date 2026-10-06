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

#include <rainbow/folded_patch_builder.hpp>

// One invocation per stored primitive; only NeedsRefinement performs quadrature.
// Slot order is irrelevant: each compact primitive receives its own lookup.
extern "C" __global__ void prepare_folded_patches(rainbow::FoldedPrepareParams p)
{
    using namespace rainbow;
    const std::uint32_t i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i>=p.patch_count)return;
    p.record_indices[i]=0xffffffffu;
    if(p.patches[i].status!=PatchCellStatus::NeedsRefinement)return;
    const auto slot=atomicAdd(p.written_count,1u);
    // Host checks written_count==capacity before publishing a successful result.
    if(slot>=p.capacity)return;
    p.record_indices[i]=slot;
    p.records[slot]=FoldedPatchBuilder::prepare(p.vertices,p.vertex_count,p.patches[i],i,
        p.incident_direction,p.with_focal?&p.focal:nullptr,p.config);
}

extern "C" __global__ void evaluate_folded_optics(rainbow::FoldedOpticsParams p)
{
    const std::uint32_t i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i>=p.wave.optical.direction_count)return;
    p.wave.optical.results[i]=rainbow::PatchOpticalEvaluator::evaluate_direction(
        p.wave.optical,i,p.with_focal?&p.wave.focal:nullptr,
        p.with_focal?p.wave.focal_results+i:nullptr,p.folded);
}

// Host load_module() requires this symbol before interpreting any v2 result or
// parameter layout. Never launched. Prevents a stale v1 fatbin from silently
// consuming the enlarged two-input structs.
extern "C" __global__ void rainbow_two_input_optics_abi_v2() {}
