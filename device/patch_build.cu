#include <rainbow/patch_geometry.hpp>
#include <cstdint>

// Existing M1 entry: preserve its name, arguments and verification pattern.
extern "C" __global__
void write_test_pattern(std::uint32_t* output,std::uint32_t count)
{
    const unsigned int index=blockIdx.x*blockDim.x+threadIdx.x;
    if(index<count) output[index]=index^0x5a5a5a5au;
}

// Pass 1: one logical cell per thread. Write every status and one small
// histogram per block. No full OutgoingPatch staging array is allocated.
extern "C" __global__
void classify_outgoing_patches(
    const rainbow::OutgoingVertex* vertices,
    rainbow::PatchBuildLayout layout,
    rainbow::PatchCellStatus* statuses,
    rainbow::PatchBlockSummary* summaries)
{
    __shared__ unsigned int histogram[rainbow::patch_status_count];
    const unsigned int lane=threadIdx.x;
    if(lane<rainbow::patch_status_count) histogram[lane]=0;
    __syncthreads();

    // Use 64 bits before the bounds check; even the final padded thread must
    // not wrap a 32-bit logical index into the beginning of the array.
    const unsigned long long index=
        static_cast<unsigned long long>(blockIdx.x)*blockDim.x+lane;
    if(index<layout.cell_count)
    {
        const auto id=static_cast<std::uint32_t>(index);
        const auto result=rainbow::PatchConstruction::make(vertices,layout,id);
        statuses[id]=result.patch.status;
        atomicAdd(&histogram[static_cast<unsigned int>(result.patch.status)],1u);
    }
    __syncthreads();
    if(lane<rainbow::patch_status_count) summaries[blockIdx.x].counts[lane]=histogram[lane];
}

// Pass 2: stable compact order = logical patch_id order.
// CPU scans block counts (not vertex data); the block-local integer scan below
// assigns deterministic ranks. Never use compact indices as stable patch IDs.
extern "C" __global__
void compact_outgoing_patches(
    const rainbow::OutgoingVertex* vertices,
    rainbow::PatchBuildLayout layout,
    const rainbow::PatchCellStatus* statuses,
    const std::uint32_t* block_offsets,
    rainbow::OutgoingPatch* patches,
    rainbow::PatchAabb* aabbs)
{
    __shared__ unsigned int prefix[rainbow::patch_threads_per_block];
    const unsigned int lane=threadIdx.x;
    const unsigned long long index=
        static_cast<unsigned long long>(blockIdx.x)*blockDim.x+lane;
    const bool active=index<layout.cell_count
        && rainbow::has_patch_geometry(statuses[static_cast<std::uint32_t>(index)]);
    prefix[lane]=active?1u:0u;
    __syncthreads();

    for(unsigned int stride=1;stride<rainbow::patch_threads_per_block;stride*=2)
    {
        const unsigned int addend=lane>=stride?prefix[lane-stride]:0u;
        __syncthreads();
        prefix[lane]+=addend;
        __syncthreads();
    }
    if(active)
    {
        const auto id=static_cast<std::uint32_t>(index);
        const std::uint32_t destination=block_offsets[blockIdx.x]+prefix[lane]-1;
        const auto result=rainbow::PatchConstruction::make(vertices,layout,id);
        patches[destination]=result.patch;
        aabbs[destination]=result.aabb;
    }
}
