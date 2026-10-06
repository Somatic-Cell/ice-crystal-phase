#include <rainbow/raindrop_tracer.hpp>
#include <rainbow/read_binary_file.hpp>
#include <rainbow/cuda_error.hpp>
#include <rainbow/optix_error.hpp>
#include <optix_stubs.h>
#include <optix_stack_size.h>
#include <algorithm>
#include <array>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <span>
#include <stdexcept>

namespace rainbow
{
namespace
{
// bounded logger．OptiX が返す必要サイズが buffer より大きい場合も読み過ぎない．
class CompileLog
{
public:
    void reset() noexcept {bytes_.fill('\0');size_=bytes_.size();}
    char* data() noexcept {return bytes_.data();}
    std::size_t* size_address() noexcept {return &size_;}
    void check(const OptixResult result,const char* operation)
    {
        const auto end=std::find(bytes_.begin(),bytes_.end(),'\0');
        const std::size_t stored=static_cast<std::size_t>(end-bytes_.begin());
        if(stored){std::fprintf(stderr,"[%s]\n",operation);std::fwrite(bytes_.data(),1,stored,stderr);std::fputc('\n',stderr);}
        if(size_>bytes_.size())std::fprintf(stderr,"[%s log truncated, required %zu bytes]\n",operation,size_);
        detail::check_optix(result,operation,__FILE__,__LINE__);
    }
private:
    std::array<char,16384> bytes_{};
    std::size_t size_=bytes_.size();
};
}
RaindropTracer::RaindropTracer(const CudaContext& cuda_context)
    :cuda_context_(cuda_context),optix_context_(cuda_context),
     raygen_record_(cuda_context),miss_record_(cuda_context),hit_record_(cuda_context),
     aabb_(cuda_context),gas_scratch_(cuda_context),gas_output_(cuda_context),
     device_params_(cuda_context),device_config_(cuda_context),vertices_(cuda_context),
     second_input_fields_(cuda_context)
{}
RaindropTracer::~RaindropTracer() noexcept {static_cast<void>(close_noexcept());}
void RaindropTracer::synchronize()
{
    RAINBOW_CUDA_CHECK(cuStreamSynchronize(cuda_context_.stream()));
    has_pending_work_=false;
}
void RaindropTracer::create_pipeline(const std::filesystem::path& optixir_path)
{
    if(is_closed_||module_) throw std::logic_error("Pipeline is already initialized or tracer was closed.");
    cuda_context_.make_current();
    const auto ir=read_binary_file(optixir_path);
    OptixModuleCompileOptions mo{};
    mo.maxRegisterCount=OPTIX_COMPILE_DEFAULT_MAX_REGISTER_COUNT;
    mo.optLevel=OPTIX_COMPILE_OPTIMIZATION_DEFAULT;
    mo.debugLevel=OPTIX_COMPILE_DEBUG_LEVEL_DEFAULT;
    OptixPipelineCompileOptions po{};
    po.traversableGraphFlags=OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_GAS;
    po.numPayloadValues=2; // local BoundaryQuery* を low/high 32bit に分ける．
    po.numAttributeValues=0;
    po.exceptionFlags=OPTIX_EXCEPTION_FLAG_NONE;
    po.pipelineLaunchParamsVariableName="raindrop_trace_params";
    po.pipelineLaunchParamsSizeInBytes=sizeof(RaindropTraceParams);
    po.usesPrimitiveTypeFlags=OPTIX_PRIMITIVE_TYPE_FLAGS_CUSTOM;
    CompileLog log;
    OptixResult result=optixModuleCreate(optix_context_.handle(),&mo,&po,ir.data(),ir.size(),
        log.data(),log.size_address(),&module_);
    log.check(result,"optixModuleCreate(raindrop)");
    OptixProgramGroupOptions go{};
    OptixProgramGroupDesc desc{};
    desc.kind=OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
    desc.raygen.module=module_;desc.raygen.entryFunctionName="__raygen__raindrop_trace_v2";
    log.reset();result=optixProgramGroupCreate(optix_context_.handle(),&desc,1,&go,log.data(),log.size_address(),&raygen_);
    log.check(result,"optixProgramGroupCreate(raygen)");
    desc={};desc.kind=OPTIX_PROGRAM_GROUP_KIND_MISS;
    desc.miss.module=module_;desc.miss.entryFunctionName="__miss__raindrop";
    log.reset();result=optixProgramGroupCreate(optix_context_.handle(),&desc,1,&go,log.data(),log.size_address(),&miss_);
    log.check(result,"optixProgramGroupCreate(miss)");
    desc={};desc.kind=OPTIX_PROGRAM_GROUP_KIND_HITGROUP;
    desc.hitgroup.moduleIS=module_;desc.hitgroup.entryFunctionNameIS="__intersection__raindrop";
    desc.hitgroup.moduleCH=module_;desc.hitgroup.entryFunctionNameCH="__closesthit__raindrop";
    log.reset();result=optixProgramGroupCreate(optix_context_.handle(),&desc,1,&go,log.data(),log.size_address(),&hitgroup_);
    log.check(result,"optixProgramGroupCreate(hitgroup)");
    OptixProgramGroup groups[]={raygen_,miss_,hitgroup_};
    OptixPipelineLinkOptions link{};
    link.maxTraceDepth=1; // 4 回の逐次 trace であり，4 段の再帰 trace ではない．
    log.reset();result=optixPipelineCreate(optix_context_.handle(),&po,&link,groups,3,log.data(),log.size_address(),&pipeline_);
    log.check(result,"optixPipelineCreate(raindrop)");
    OptixStackSizes stacks{};
    for(auto group:groups) RAINBOW_OPTIX_CHECK(optixUtilAccumulateStackSizes(group,&stacks,pipeline_));
    unsigned traversal=0,state=0,continuation=0;
    RAINBOW_OPTIX_CHECK(optixUtilComputeStackSizes(&stacks,1,0,0,&traversal,&state,&continuation));
    RAINBOW_OPTIX_CHECK(optixPipelineSetStackSize(pipeline_,traversal,state,continuation,1));
    RAINBOW_OPTIX_CHECK(optixSbtRecordPackHeader(raygen_,&host_raygen_));
    RAINBOW_OPTIX_CHECK(optixSbtRecordPackHeader(miss_,&host_miss_));
    RAINBOW_OPTIX_CHECK(optixSbtRecordPackHeader(hitgroup_,&host_hit_));
    raygen_record_.allocate(1);miss_record_.allocate(1);hit_record_.allocate(1);
    aabb_.allocate(1);device_params_.allocate(1);device_config_.allocate(1);
    has_pending_work_=true;
    const auto stream=cuda_context_.stream();
    raygen_record_.upload_async(std::span<const HeaderRecord>{&host_raygen_,1},stream);
    miss_record_.upload_async(std::span<const HeaderRecord>{&host_miss_,1},stream);
    hit_record_.upload_async(std::span<const HeaderRecord>{&host_hit_,1},stream);
    synchronize();
}
void RaindropTracer::build_drop_gas()
{
    // normalized drop 一個．shape の変更後でも bound を更新し，同じ型の GAS を作る．
    const float r=config_.shape.outer_radius+0x1p-15f;
    host_aabb_={-r,-r,-r,r,r,r};
    const auto stream=cuda_context_.stream();
    has_pending_work_=true;
    aabb_.upload_async(std::span<const OptixAabb>{&host_aabb_,1},stream);
    CUdeviceptr address=aabb_.address();
    unsigned flags=OPTIX_GEOMETRY_FLAG_DISABLE_ANYHIT;
    OptixBuildInput input{};
    input.type=OPTIX_BUILD_INPUT_TYPE_CUSTOM_PRIMITIVES;
    input.customPrimitiveArray.aabbBuffers=&address;
    input.customPrimitiveArray.numPrimitives=1;
    input.customPrimitiveArray.strideInBytes=sizeof(OptixAabb);
    input.customPrimitiveArray.flags=&flags;
    input.customPrimitiveArray.numSbtRecords=1;
    OptixAccelBuildOptions options{};
    options.buildFlags=OPTIX_BUILD_FLAG_PREFER_FAST_TRACE;
    options.operation=OPTIX_BUILD_OPERATION_BUILD;
    OptixAccelBufferSizes sizes{};
    RAINBOW_OPTIX_CHECK(optixAccelComputeMemoryUsage(optix_context_.handle(),&options,&input,1,&sizes));
    // trace() 冒頭で前回 launch は同期済み．古い GAS を再利用中に解放しない．
    gas_scratch_.close();gas_output_.close();
    gas_scratch_.allocate(sizes.tempSizeInBytes);gas_output_.allocate(sizes.outputSizeInBytes);
    OptixTraversableHandle handle=0;
    RAINBOW_OPTIX_CHECK(optixAccelBuild(optix_context_.handle(),stream,&options,&input,1,
        gas_scratch_.address(),gas_scratch_.byte_size(),gas_output_.address(),gas_output_.byte_size(),&handle,nullptr,0));
    params_.traversable=static_cast<std::uint64_t>(handle);
    // build と launch は同じ stream．ここに追加の同期は不要．
}
void RaindropTracer::trace(const RaindropSettings& settings) { trace_impl(settings, false); }
void RaindropTracer::trace_unpolarized(const RaindropSettings& settings) { trace_impl(settings, true); }
void RaindropTracer::trace_impl(const RaindropSettings& settings, const bool unpolarized)
{
    if(is_closed_||!pipeline_) throw std::logic_error("Call create_pipeline before trace.");
    cuda_context_.make_current();
    if(has_pending_work_) synchronize();
    has_result_=false;
    auto effective_settings = settings;
    if(unpolarized) effective_settings.incident_field = {{1,0},{0,0}}; // unit column-0 descriptor only
    config_=effective_settings.make_config();
    unpolarized_=unpolarized;
    const std::size_t total=std::size_t(config_.vertex_count())*4;
    if(vertices_.element_count()!=total){vertices_.close();vertices_.allocate(total);}
    if(unpolarized_)
    {
        if(second_input_fields_.element_count()!=total)
        { second_input_fields_.close(); second_input_fields_.allocate(total); }
    }
    else second_input_fields_.close();
    params_.second_input_fields=unpolarized_ ? second_input_fields_.data() : nullptr;
    params_.vertices=vertices_.data();
    params_.config=device_config_.data();
    build_drop_gas();
    const auto stream=cuda_context_.stream();
    device_config_.upload_async(std::span<const RaindropTraceConfig>{&config_,1},stream);
    device_params_.upload_async(std::span<const RaindropTraceParams>{&params_,1},stream);
    OptixShaderBindingTable sbt{};
    sbt.raygenRecord=raygen_record_.address();
    sbt.missRecordBase=miss_record_.address();sbt.missRecordCount=1;sbt.missRecordStrideInBytes=sizeof(HeaderRecord);
    sbt.hitgroupRecordBase=hit_record_.address();sbt.hitgroupRecordCount=1;sbt.hitgroupRecordStrideInBytes=sizeof(HeaderRecord);
    RAINBOW_OPTIX_CHECK(optixLaunch(pipeline_,stream,device_params_.address(),sizeof(params_),&sbt,
        config_.vertex_count(),1,1));
    synchronize();
    has_result_=true;
}
std::vector<OutgoingVertex> RaindropTracer::download_vertices() const
{
    if(!has_result_) throw std::logic_error("No completed trace result.");
    std::vector<OutgoingVertex> output(vertices_.element_count());
    vertices_.download(std::span<OutgoingVertex>{output});
    return output;
}
std::vector<Field32> RaindropTracer::download_second_input_fields() const
{
    if(!has_result_ || !unpolarized_) throw std::logic_error("No completed two-input trace.");
    std::vector<Field32> output(second_input_fields_.element_count());
    second_input_fields_.download(std::span<Field32>{output});
    return output;
}
void RaindropTracer::write_csv(const std::filesystem::path& output_path) const
{
    const auto output=download_vertices();
    write_csv(output_path,std::span<const OutgoingVertex>{output});
}
void RaindropTracer::write_csv(const std::filesystem::path& output_path,
    const std::span<const OutgoingVertex> output) const
{
    if(!has_result_ || output.size()!=vertices_.element_count())
        throw std::invalid_argument("CSV data must match the completed trace.");
    const auto second = unpolarized_ ? download_second_input_fields() : std::vector<Field32>{};
    std::ofstream file(output_path);
    if(!file) throw std::runtime_error("Cannot open outgoing vertex CSV.");
    file<<std::setprecision(std::numeric_limits<float>::max_digits10);
    const auto& c=config_;
    file<<"# format="<<(unpolarized_?"rainbow_outgoing_vertices_v2":"rainbow_outgoing_vertices_v1")<<"\n# radius_mm="<<c.radius_mm<<"\n# wavelength_nm="<<c.wavelength_nm
        <<"\n# interior_index="<<c.interior_index<<"\n# exterior_index="<<c.exterior_index
        <<"\n# grid="<<c.grid_width<<","<<c.grid_height<<"\n# grid_half_extent_drop="<<c.grid_half_extent
        <<"\n# reference_distances_drop="<<c.reference_distance<<","<<c.outgoing_reference_distance
        <<"\n# incident_direction="<<c.incident_direction.x<<","<<c.incident_direction.y<<","<<c.incident_direction.z
        <<"\n# incident_basis_x="<<c.incident_basis_x.x<<","<<c.incident_basis_x.y<<","<<c.incident_basis_x.z;
    if(unpolarized_)
        file << "\n# input_states=unpolarized_two_orthogonal_unit_Jones_inputs"
             << "\n# incident_polarization=unpolarized\n# incident_total_intensity=1"
             << "\n# input_coherency=0.5,0,0,0.5"
             << "\n# response_column_0_input=1,0,0,0\n# response_column_1_input=0,0,1,0"
             << "\n# field_columns=two_unit_input_responses_NOT_one_unpolarized_field";
    else file << "\n# incident_field="<<c.incident_field.x.real<<","<<c.incident_field.x.imag
              <<","<<c.incident_field.y.real<<","<<c.incident_field.y.imag;
    file << "\n# coefficients=";
    for(unsigned k=0;k<8;++k) file<<(k?",":"")<<c.shape.coefficients[k];
    file<<"\n# coordinates: +y up; polar theta from -y; position normalized by radius.\n# field excludes propagation/focal/patch-area factors.\n"
        <<"family,grid_index,ix,iy,status,diagnostics,px_drop,py_drop,pz_drop,wx,wy,wz,bx,by,bz,";
    if(unpolarized_) file<<"J00_real,J00_imag,J10_real,J10_imag,turns,fraction,J01_real,J01_imag,J11_real,J11_imag\n";
    else file<<"Ex_real,Ex_imag,Ey_real,Ey_imag,turns,fraction\n";
    constexpr const char* families[]={"R","TT","TRT","TRRT"};
    for(unsigned p=0;p<4;++p) for(std::uint32_t i=0;i<c.vertex_count();++i)
    {
        const auto& v=output[std::size_t(p)*c.vertex_count()+i];
        file<<families[p]<<','<<i<<','<<i%c.grid_width<<','<<i/c.grid_width<<','<<static_cast<unsigned>(v.status)<<','<<v.diagnostics<<','
            <<v.position_drop.x<<','<<v.position_drop.y<<','<<v.position_drop.z<<','
            <<v.direction_drop.x<<','<<v.direction_drop.y<<','<<v.direction_drop.z<<','
            <<v.basis_x.x<<','<<v.basis_x.y<<','<<v.basis_x.z<<','
            <<v.field.x.real<<','<<v.field.x.imag<<','<<v.field.y.real<<','<<v.field.y.imag<<','
            <<v.optical_cycles.turns<<','<<v.optical_cycles.fraction;
        if(unpolarized_)
        {
            const auto& f=second[std::size_t(p)*c.vertex_count()+i];
            file<<','<<f.x.real<<','<<f.x.imag<<','<<f.y.real<<','<<f.y.imag;
        }
        file<<'\n';
    }
    file.flush();if(!file)throw std::runtime_error("Writing outgoing vertex CSV failed.");
}
bool RaindropTracer::close_noexcept() noexcept
{
    if(is_closed_)return true;
    if(!detail::report_cuda_cleanup_result(cuCtxSetCurrent(cuda_context_.handle()),"cuCtxSetCurrent(tracer cleanup)"))return false;
    bool ok=true;
    if(has_pending_work_){ok=detail::report_cuda_cleanup_result(cuStreamSynchronize(cuda_context_.stream()),"cuStreamSynchronize(tracer cleanup)")&&ok;has_pending_work_=false;}
    ok=second_input_fields_.close_noexcept()&&ok;
    ok=vertices_.close_noexcept()&&ok;ok=device_params_.close_noexcept()&&ok;
    ok=device_config_.close_noexcept()&&ok;
    ok=gas_scratch_.close_noexcept()&&ok;ok=gas_output_.close_noexcept()&&ok;ok=aabb_.close_noexcept()&&ok;
    ok=raygen_record_.close_noexcept()&&ok;ok=miss_record_.close_noexcept()&&ok;ok=hit_record_.close_noexcept()&&ok;
    if(pipeline_){ok=detail::report_optix_cleanup_result(optixPipelineDestroy(pipeline_),"optixPipelineDestroy")&&ok;pipeline_=nullptr;}
    for(auto* group:{&hitgroup_,&miss_,&raygen_}) if(*group){ok=detail::report_optix_cleanup_result(optixProgramGroupDestroy(*group),"optixProgramGroupDestroy")&&ok;*group=nullptr;}
    if(module_){ok=detail::report_optix_cleanup_result(optixModuleDestroy(module_),"optixModuleDestroy")&&ok;module_=nullptr;}
    ok=optix_context_.close_noexcept()&&ok;
    is_closed_=true;has_result_=false;
    return ok;
}
void RaindropTracer::close(){if(!close_noexcept())throw std::runtime_error("Raindrop tracer cleanup failed; see stderr.");}
}
