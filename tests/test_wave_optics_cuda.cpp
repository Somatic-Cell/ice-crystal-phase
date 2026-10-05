#include "wave_optics_test_cases.hpp"
#include <rainbow/patch_optics.hpp>
#include <rainbow/patch_query.hpp>
#include <rainbow/cuda_error.hpp>
#include <rainbow/cuda_module.hpp>
#include <rainbow/device_buffer.hpp>
#include <cstdlib>
#include <filesystem>
#include <iostream>

namespace
{
using namespace rainbow;using namespace rainbow::tests;
// Test-only guard, declared after device buffers so pending launches finish
// before stack unwinding frees input/output storage. Not a runtime substitute.
class SynchronizeOnExit
{
public:
    explicit SynchronizeOnExit(const CudaContext& c):c_(c){}
    ~SynchronizeOnExit() noexcept
    {
        static_cast<void>(detail::report_cuda_cleanup_result(cuCtxSetCurrent(c_.handle()),"wave-test context"));
        static_cast<void>(detail::report_cuda_cleanup_result(cuStreamSynchronize(c_.stream()),"wave-test sync"));
    }
    SynchronizeOnExit(const SynchronizeOnExit&)=delete;
    SynchronizeOnExit& operator=(const SynchronizeOnExit&)=delete;
private:const CudaContext& c_;
};
template<class T> void upload(DeviceBuffer<T>& buffer,const std::vector<T>& data,const CudaContext& cuda)
{
    buffer.allocate(data.size());buffer.upload_async(std::span<const T>{data},cuda.stream());
    RAINBOW_CUDA_CHECK(cuStreamSynchronize(cuda.stream()));
}
template<class P> void launch(CudaModule& module,const char* name,P& params,std::uint32_t count,const CudaContext& cuda)
{
    void* args[]={&params};constexpr unsigned threads=128;
    RAINBOW_CUDA_CHECK(cuLaunchKernel(module.find_function(name),count/threads+(count%threads!=0u?1u:0u),1,1,
        threads,1,1,0,cuda.stream(),args,nullptr));
}
void numerical_tests(const CudaContext& cuda,const std::filesystem::path& module_path)
{
    CudaModule module(cuda);module.load_fatbin(module_path);
    auto f=make_wave_fixture();auto& h=f.optical;
    DeviceBuffer<OutgoingVertex> vertices(cuda);DeviceBuffer<OutgoingPatch> patches(cuda);
    DeviceBuffer<Vec3> directions(cuda);DeviceBuffer<std::uint64_t> offsets(cuda);
    DeviceBuffer<PatchQuerySummary> summaries(cuda);DeviceBuffer<PatchQueryHit> hits(cuda);
    DeviceBuffer<PatchOpticalResult> raw(cuda);DeviceBuffer<FocalOpticalResult> focal(cuda);
    std::vector<PatchOpticalResult> host_raw(h.directions.size()+2);
    std::vector<FocalOpticalResult> host_focal(h.directions.size()+2);
    host_raw.front().reserved=0x12345678u;host_raw.back().reserved=0x87654321u;
    host_focal.front().flags=0x12345678u;host_focal.back().flags=0x87654321u;
    SynchronizeOnExit wait(cuda);
    upload(vertices,h.vertices,cuda);upload(patches,h.patches,cuda);upload(directions,h.directions,cuda);
    upload(offsets,h.offsets,cuda);upload(summaries,h.summaries,cuda);upload(hits,h.hits,cuda);
    upload(raw,host_raw,cuda);upload(focal,host_focal,cuda);
    auto p=h.view();p.vertices=vertices.data();p.patches=patches.data();p.directions=directions.data();
    p.offsets=offsets.data();p.summaries=summaries.data();p.hits=hits.data();p.results=raw.data()+1;
    WaveOpticsParams w{p,f.focal,focal.data()+1};
    launch(module,"evaluate_focal_optics",w,p.direction_count,cuda);
    RAINBOW_CUDA_CHECK(cuStreamSynchronize(cuda.stream()));raw.download(host_raw);focal.download(host_focal);
    require_optics(host_raw.front().reserved==0x12345678u&&host_raw.back().reserved==0x87654321u
        &&host_focal.front().flags==0x12345678u&&host_focal.back().flags==0x87654321u,"Wave CUDA guard overwritten.");
    for(std::uint32_t i=0;i<p.direction_count;++i)
    {
        FocalOpticalResult reference{};
        const auto original=PatchOpticalEvaluator::evaluate_direction(h.view(),i,&f.focal,&reference);
        compare_optics_results(host_raw[i+1],original);compare_focal(host_focal[i+1],reference);
    }
    std::cout<<"CUDA focal evaluator: "<<p.direction_count<<" directions passed.\n";
    auto grid=make_diffraction_fixture(120,80);
    DeviceBuffer<FocalOpticalResult> source(cuda);DeviceBuffer<RainbowTransition> edges(cuda);
    DeviceBuffer<DiffractionResult> output(cuda);
    std::vector<RainbowTransition> host_edges(grid.input.size()+2);
    std::vector<DiffractionResult> host_out(grid.input.size()+2);
    host_edges.front().flags=0x12345678u;host_edges.back().flags=0x87654321u;
    host_out.front().flags=0x12345678u;host_out.back().flags=0x87654321u;
    SynchronizeOnExit wait_filter(cuda);
    upload(source,grid.input,cuda);upload(edges,host_edges,cuda);upload(output,host_out,cuda);
    DiffractionParams d{source.data(),edges.data()+1,output.data()+1,grid.rows,grid.columns,grid.config};
    const auto n=grid.rows*grid.columns;
    launch(module,"detect_rainbow_transitions",d,n,cuda);
    launch(module,"filter_rainbow_diffraction",d,n,cuda);
    RAINBOW_CUDA_CHECK(cuStreamSynchronize(cuda.stream()));edges.download(host_edges);output.download(host_out);
    require_optics(host_edges.front().flags==0x12345678u&&host_edges.back().flags==0x87654321u
        &&host_out.front().flags==0x12345678u&&host_out.back().flags==0x87654321u,"Diffraction CUDA guard overwritten.");
    auto ref=grid.view();
    for(std::uint32_t i=0;i<n;++i)
    {
        require_optics(host_edges[i+1].flags==grid.edges[i].flags,"CUDA edge detection mismatch.");
        compare_diffraction(host_out[i+1],RainbowDiffraction::filter(ref,i));
    }
    std::cout<<"CUDA diffraction: "<<n<<" directions passed.\n";
}
void integration_test(const CudaContext& cuda,const std::filesystem::path& optical,
                      const std::filesystem::path& build,const std::filesystem::path& query_path)
{
    // Small, genuinely two-dimensional chord patches, not a duplicate ray solver.
    RaindropTraceConfig config{};config.grid_width=config.grid_height=5;config.grid_half_extent=1;
    config.incident_direction={0,0,1};config.incident_basis_x={1,0,0};config.radius_mm=.4f;
    std::vector<OutgoingVertex> host(100);
    for(unsigned f=0;f<4;++f)for(unsigned y=0;y<5;++y)for(unsigned x=0;x<5;++x)
    {
        auto& v=host[f*25+y*5+x];v.status=VertexStatus::Valid;
        v.direction_drop=Vec3{.1f*(float(x)-2),.1f*(float(y)-2),1}.normalized();
        v.basis_x=TransverseFrame::from_direction(v.direction_drop).e0;
        v.field={{1,0},{0,0}};v.optical_cycles={1000,0};
    }
    OptixContext optix(cuda);DeviceBuffer<OutgoingVertex> vertices(cuda);upload(vertices,host,cuda);
    PatchAccel accel(cuda,optix);accel.load_module(build);accel.build(vertices,config);
    PatchQuery query(cuda,optix);query.create_pipeline(query_path);query.query_grid(accel,config,18,36);
    PatchOptics optics(cuda);optics.load_module(optical);
    WaveOpticsSettings settings{};optics.evaluate_wave(accel,query,config,settings);
    require_optics(optics.has_wave_result()&&optics.wave_statistics().focal_errors==0,"Wave production API failed.");
    require_optics(optics.wave_statistics().corrected_hits>0,"Integration test had no optical hits.");
    const auto first=std::vector<FocalOpticalResult>(optics.host_focal_results().begin(),optics.host_focal_results().end());
    const auto original=std::vector<PatchOpticalResult>(optics.host_results().begin(),optics.host_results().end());
    // Every family in this test has the same positive derivative: common +pi/2.
    for(std::size_t i=0;i<first.size();++i)if(first[i].valid())
        require_optics(wave_near(first[i].intensity_s,original[i].regular_partial_path_s)
            &&wave_near(first[i].intensity_p,original[i].regular_partial_path_p),"Common correction changed intensity.");
    optics.evaluate_wave(accel,query,config,settings);
    for(std::size_t i=0;i<first.size();++i)compare_focal(first[i],optics.host_focal_results()[i]);
    optics.evaluate(accel,query,config);
    require_optics(!optics.has_wave_result(),"Stale wave output after path-only evaluate.");
    for(std::size_t i=0;i<original.size();++i)compare_optics_results(original[i],optics.host_results()[i]);
    optics.close();query.close();accel.close();vertices.close();
    std::cout<<"Wave production integration: GAS/query -> focal/filter; repeat and legacy reuse passed.\n";
}
}
#if defined(_WIN32)
int wmain(int argc,wchar_t* argv[])
#else
int main(int argc,char* argv[])
#endif
{
    try
    {
        if(argc!=4)throw std::invalid_argument("Usage: wave_test <patch_optics.fatbin> <patch_build.fatbin> <patch_query.optixir>");
        CudaContext cuda{0};numerical_tests(cuda,std::filesystem::path(argv[1]));
        integration_test(cuda,std::filesystem::path(argv[1]),std::filesystem::path(argv[2]),std::filesystem::path(argv[3]));
        return EXIT_SUCCESS;
    }
    catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<'\n';return EXIT_FAILURE;}
}
