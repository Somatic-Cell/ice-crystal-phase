#include <rainbow/patch_accel.hpp>
#include <rainbow/raindrop_tracer.hpp>
#include <rainbow/cuda_error.hpp>
#include <rainbow/optix_error.hpp>
#include <optix_stubs.h>

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <span>
#include <stdexcept>

namespace rainbow
{
// OptiX consumes the same 6 floats the CUDA kernel wrote. Verify all offsets,
// not only sizeof. DeviceBuffer's cuMemAlloc satisfies the address alignment.
static_assert(sizeof(PatchAabb)==sizeof(OptixAabb));
static_assert(offsetof(PatchAabb,min_x)==offsetof(OptixAabb,minX));
static_assert(offsetof(PatchAabb,min_y)==offsetof(OptixAabb,minY));
static_assert(offsetof(PatchAabb,min_z)==offsetof(OptixAabb,minZ));
static_assert(offsetof(PatchAabb,max_x)==offsetof(OptixAabb,maxX));
static_assert(offsetof(PatchAabb,max_y)==offsetof(OptixAabb,maxY));
static_assert(offsetof(PatchAabb,max_z)==offsetof(OptixAabb,maxZ));

PatchAccel::PatchAccel(const CudaContext& cuda_context,const OptixContext& optix_context) noexcept
    :cuda_context_(cuda_context),optix_context_(optix_context),module_(cuda_context),
     statuses_(cuda_context),block_summaries_(cuda_context),block_offsets_(cuda_context),
     patches_(cuda_context),aabbs_(cuda_context),gas_scratch_(cuda_context),gas_output_(cuda_context)
{}
PatchAccel::~PatchAccel() noexcept {static_cast<void>(close_noexcept());}

void PatchAccel::load_module(const std::filesystem::path& fatbin_path)
{
    if(is_closed_) throw std::logic_error("PatchAccel is closed.");
    module_.load_fatbin(fatbin_path);
}
void PatchAccel::synchronize()
{
    RAINBOW_CUDA_CHECK(cuStreamSynchronize(cuda_context_.stream()));
    has_pending_work_=false;
}
void PatchAccel::clear_result()
{
    has_result_=false;
    traversable_=0;
    source_vertices_=0;
    source_second_fields_=0;
    statistics_={};
    // Called only after completion of the owning stream. No per-buffer hidden
    // synchronizations are added to DeviceBuffer.
    gas_output_.close(); gas_scratch_.close();
    aabbs_.close(); patches_.close();
    statuses_.close(); block_offsets_.close(); block_summaries_.close();
}

void PatchAccel::build(const RaindropTracer& tracer)
{
    if(!tracer.has_result()) throw std::invalid_argument("Patch build requires a completed trace.");
    build(tracer.vertices(), tracer.config(), tracer.is_unpolarized() ? &tracer.second_input_fields() : nullptr);
}
void PatchAccel::build(const DeviceBuffer<OutgoingVertex>& vertices,const RaindropTraceConfig& config,
                      const DeviceBuffer<Field32>* second_input_fields)
{
    if(is_closed_ || !module_.is_loaded()) throw std::logic_error("Load the patch module before build().");
    if(optix_context_.handle()==nullptr) throw std::logic_error("OptiX context is not live.");
    PatchBuildLayout next_layout{};
    if(!PatchBuildLayout::try_make(config,next_layout))
        throw std::invalid_argument("Invalid patch grid or cell area.");
    if(vertices.element_count()!=std::size_t(next_layout.vertices_per_path)*4)
        throw std::invalid_argument("Patch source must contain all four path-major vertex grids.");
    if(second_input_fields && (second_input_fields->element_count()!=vertices.element_count()
                              || second_input_fields->address()==0))
        throw std::invalid_argument("Unpolarized source needs a complete second response column.");

    cuda_context_.make_current();
    synchronize();
    clear_result();
    layout_=next_layout;
    source_vertices_=vertices.address();
    source_second_fields_=second_input_fields ? second_input_fields->address() : 0;
    const unsigned int block_count=layout_.cell_count/patch_threads_per_block
        +(layout_.cell_count%patch_threads_per_block!=0?1u:0u);
    statuses_.allocate(layout_.cell_count);
    block_summaries_.allocate(block_count);
    block_offsets_.allocate(block_count);
    host_summaries_.resize(block_count);
    host_offsets_.resize(block_count);

    const CUstream stream=cuda_context_.stream();
    CUdeviceptr source=source_vertices_,statuses=statuses_.address(),summaries=block_summaries_.address();
    void* classify_arguments[]={&source,&layout_,&statuses,&summaries};
    const CUfunction classify=module_.find_function("classify_outgoing_patches");
    const CUfunction compact=module_.find_function("compact_outgoing_patches");
    has_pending_work_=true;
    RAINBOW_CUDA_CHECK(cuLaunchKernel(classify,block_count,1,1,patch_threads_per_block,1,1,
        0,stream,classify_arguments,nullptr));
    synchronize();
    block_summaries_.download(std::span<PatchBlockSummary>{host_summaries_});

    // Only O(cell_count/256) integers cross to the host here. The exclusive
    // scan is exact and deterministic; it does not accumulate optical energy.
    std::uint64_t active_count=0,total=0;
    for(unsigned int block=0;block<block_count;++block)
    {
        host_offsets_[block]=static_cast<std::uint32_t>(active_count);
        const auto& summary=host_summaries_[block];
        for(unsigned int s=0;s<patch_status_count;++s)
        {
            statistics_.cells[s]+=summary.counts[s];
            total+=summary.counts[s];
            if(has_patch_geometry(static_cast<PatchCellStatus>(s))) active_count+=summary.counts[s];
        }
        if(active_count>layout_.cell_count)
            throw std::runtime_error("Patch block count is corrupt.");
    }
    statistics_.logical_cell_count=layout_.cell_count;
    if(total!=layout_.cell_count || statistics_.count(PatchCellStatus::Unwritten)!=0)
        throw std::runtime_error("Not all logical patch cells were classified.");
    if(statistics_.count(PatchCellStatus::InvalidVertex)!=0
       || statistics_.count(PatchCellStatus::InvalidGeometry)!=0)
        throw std::runtime_error("Invalid source vertices: patch GAS was NOT built. Cell statuses remain available.");

    statistics_.patch_count=static_cast<std::uint32_t>(active_count);
    if(active_count==0)
    {
        // Empty GAS: handle 0 is intentional. A future query must skip optixTrace.
        has_result_=true;
        return;
    }
    patches_.allocate(static_cast<std::size_t>(active_count));
    aabbs_.allocate(static_cast<std::size_t>(active_count));
    has_pending_work_=true;
    block_offsets_.upload_async(std::span<const std::uint32_t>{host_offsets_},stream);
    CUdeviceptr offsets=block_offsets_.address(),patches=patches_.address(),bounds=aabbs_.address();
    void* compact_arguments[]={&source,&layout_,&statuses,&offsets,&patches,&bounds};
    RAINBOW_CUDA_CHECK(cuLaunchKernel(compact,block_count,1,1,patch_threads_per_block,1,1,
        0,stream,compact_arguments,nullptr));
    // GAS build follows compaction on exactly the same stream: no readback of AABBs.
    build_gas();
    synchronize();
    // Scratch is not part of the traversable's lifetime. Release it after build.
    gas_scratch_.close();
    has_result_=true;
}

void PatchAccel::build_gas()
{
    unsigned int maximum_primitives=0;
    RAINBOW_OPTIX_CHECK(optixDeviceContextGetProperty(optix_context_.handle(),
        OPTIX_DEVICE_PROPERTY_LIMIT_MAX_PRIMITIVES_PER_GAS,&maximum_primitives,sizeof(maximum_primitives)));
    if(statistics_.patch_count>maximum_primitives)
        throw std::length_error("Patch count exceeds the device's single-GAS primitive limit.");
    if(aabbs_.address()%OPTIX_AABB_BUFFER_BYTE_ALIGNMENT!=0)
        throw std::runtime_error("Patch AABB device address is not aligned.");
    CUdeviceptr address=aabbs_.address();
    // Keep any-hit enabled. In the next stage, all intersections must be
    // enumerated; spatial splitting must not duplicate an any-hit invocation.
    unsigned int flags=OPTIX_GEOMETRY_FLAG_REQUIRE_SINGLE_ANYHIT_CALL;
    OptixBuildInput input{};
    input.type=OPTIX_BUILD_INPUT_TYPE_CUSTOM_PRIMITIVES;
    input.customPrimitiveArray.aabbBuffers=&address;
    input.customPrimitiveArray.numPrimitives=statistics_.patch_count;
    input.customPrimitiveArray.strideInBytes=sizeof(PatchAabb);
    input.customPrimitiveArray.flags=&flags;
    input.customPrimitiveArray.numSbtRecords=1;

    OptixAccelBuildOptions options{};
    options.buildFlags=OPTIX_BUILD_FLAG_PREFER_FAST_TRACE;
    options.operation=OPTIX_BUILD_OPERATION_BUILD;
    OptixAccelBufferSizes sizes{};
    RAINBOW_OPTIX_CHECK(optixAccelComputeMemoryUsage(optix_context_.handle(),&options,&input,1,&sizes));
    gas_scratch_.allocate(sizes.tempSizeInBytes);
    gas_output_.allocate(sizes.outputSizeInBytes);
    statistics_.gas_byte_size=sizes.outputSizeInBytes;
    RAINBOW_OPTIX_CHECK(optixAccelBuild(optix_context_.handle(),cuda_context_.stream(),&options,&input,1,
        gas_scratch_.address(),gas_scratch_.byte_size(),gas_output_.address(),gas_output_.byte_size(),
        &traversable_,nullptr,0));
}

void PatchAccel::write_csv(const std::filesystem::path& path) const
{
    if(!has_result_) throw std::logic_error("No completed patch build.");
    std::vector<OutgoingPatch> records(patches_.element_count());
    std::vector<PatchAabb> bounds(aabbs_.element_count());
    std::vector<PatchCellStatus> cells(statuses_.element_count());
    patches_.download(std::span<OutgoingPatch>{records});
    aabbs_.download(std::span<PatchAabb>{bounds});
    statuses_.download(std::span<PatchCellStatus>{cells});
    std::ofstream output(path);
    if(!output) throw std::runtime_error("Cannot create patch CSV.");
    output<<std::setprecision(std::numeric_limits<float>::max_digits10);
    output<<"# format=rainbow_patch_geometry_v1\n"
          <<"# NOT a phase-function LUT; boundary completion and fold refinement are pending.\n"
          <<"# source corner order=00,10,01,11; geometry=P(u,v), NOT normalized P.\n"
          <<"# area unit=radius^2; physical area_m2=incident_area_drop2*(radius_mm*1e-3)^2.\n"
          <<"patch_id,family,cell_x,cell_y,status,compact_index,v00,v10,v01,v11,"
          <<"incident_area_drop2,signed_solid_angle_sr,min_x,min_y,min_z,max_x,max_y,max_z\n";
    std::size_t compact_index=0;
    for(std::uint32_t id=0;id<layout_.cell_count;++id)
    {
        std::uint32_t ids[4]; layout_.corner_indices(id,ids);
        const auto status=cells[id];
        const std::uint32_t cell=id%layout_.cells_per_path;
        output<<id<<','<<id/layout_.cells_per_path<<','<<cell%(layout_.grid_width-1)
              <<','<<cell/(layout_.grid_width-1)<<','<<static_cast<unsigned>(status)<<',';
        if(has_patch_geometry(status)) output<<compact_index;
        else output<<-1;
        for(const auto vertex:ids) output<<','<<vertex;
        output<<','<<layout_.incident_area_drop2;
        if(has_patch_geometry(status))
        {
            if(compact_index>=records.size() || records[compact_index].patch_id!=id
               || records[compact_index].status!=status)
                throw std::runtime_error("Patch compaction mapping is inconsistent.");
            const auto& b=bounds[compact_index];
            output<<','<<records[compact_index].signed_solid_angle_sr
                  <<','<<b.min_x<<','<<b.min_y<<','<<b.min_z
                  <<','<<b.max_x<<','<<b.max_y<<','<<b.max_z;
            ++compact_index;
        }
        else output<<",nan,nan,nan,nan,nan,nan,nan";
        output<<'\n';
    }
    if(compact_index!=records.size()) throw std::runtime_error("Patch count is inconsistent.");
    output.flush();
    if(!output) throw std::runtime_error("Writing patch CSV failed.");
}

bool PatchAccel::close_noexcept() noexcept
{
    if(is_closed_) return true;
    if(!detail::report_cuda_cleanup_result(cuCtxSetCurrent(cuda_context_.handle()),"cuCtxSetCurrent(patch cleanup)"))
        return false;
    bool succeeded=true;
    if(has_pending_work_)
    {
        succeeded=detail::report_cuda_cleanup_result(cuStreamSynchronize(cuda_context_.stream()),"cuStreamSynchronize(patch cleanup)");
        has_pending_work_=false;
    }
    succeeded=gas_output_.close_noexcept()&&succeeded;
    succeeded=gas_scratch_.close_noexcept()&&succeeded;
    succeeded=aabbs_.close_noexcept()&&succeeded;
    succeeded=patches_.close_noexcept()&&succeeded;
    succeeded=statuses_.close_noexcept()&&succeeded;
    succeeded=block_offsets_.close_noexcept()&&succeeded;
    succeeded=block_summaries_.close_noexcept()&&succeeded;
    succeeded=module_.close_noexcept()&&succeeded;
    traversable_=0;source_vertices_=0;has_result_=false;is_closed_=true;
    return succeeded;
}
void PatchAccel::close()
{
    if(!close_noexcept()) throw std::runtime_error("PatchAccel cleanup failed. See stderr.");
}

} // namespace rainbow
