#include "folded_test_support.hpp"
#include <rainbow/patch_optics.hpp>
#include <rainbow/patch_query.hpp>
#include <rainbow/cuda_error.hpp>
#include <rainbow/cuda_module.hpp>
#include <rainbow/device_buffer.hpp>
#include <filesystem>
#include <iostream>

namespace
{
using namespace rainbow;using namespace folded_tests;
struct Wait
{
    const CudaContext& cuda;
    ~Wait(){static_cast<void>(detail::report_cuda_cleanup_result(cuCtxSetCurrent(cuda.handle()),"fold test context"));
        static_cast<void>(detail::report_cuda_cleanup_result(cuStreamSynchronize(cuda.stream()),"fold test sync"));}
};
template<class T> void upload(DeviceBuffer<T>& b,const std::vector<T>& a,const CudaContext& c)
{b.allocate(a.size());b.upload_async(std::span<const T>{a},c.stream());RAINBOW_CUDA_CHECK(cuStreamSynchronize(c.stream()));}
template<class P> void launch(CudaModule& m,const char* name,P& p,unsigned n,const CudaContext& c,unsigned threads=32)
{
    void* args[]={&p};
    RAINBOW_CUDA_CHECK(cuLaunchKernel(m.find_function(name),n/threads+(n%threads!=0u),1,1,threads,1,1,0,c.stream(),args,nullptr));
}
void witness(const CudaContext& c,const std::filesystem::path& module_path)
{
    auto a=make_witness();auto host=prepare(a);CudaModule module(c);module.load_fatbin(module_path);
    DeviceBuffer<OutgoingVertex> v(c);DeviceBuffer<OutgoingPatch> patches(c);
    DeviceBuffer<Vec3> directions(c);DeviceBuffer<PatchQueryHit> hits(c);
    DeviceBuffer<PatchQuerySummary> summaries(c);DeviceBuffer<std::uint64_t> offsets(c);
    DeviceBuffer<std::uint32_t> lookup(c),count(c);DeviceBuffer<FoldedPatchRecord> records(c);
    DeviceBuffer<PatchOpticalResult> result(c);DeviceBuffer<FocalOpticalResult> focal(c);
    const auto n=static_cast<unsigned>(a.directions.size());const auto capacity=static_cast<unsigned>(host.records.size());
    std::vector<std::uint32_t> hm(a.patches.size()+2,0xabcdef01u),hc={0xabcdef01u,0,0xabcdef01u};
    std::vector<FoldedPatchRecord> hr(capacity+2);hr.front().flags=hr.back().flags=0xabcdef01u;
    std::vector<PatchOpticalResult> ho(n+2);ho.front().reserved=ho.back().reserved=0xabcdef01u;
    std::vector<FocalOpticalResult> hf(n+2);hf.front().flags=hf.back().flags=0xabcdef01u;
    Wait wait{c};
    upload(v,a.vertices,c);upload(patches,a.patches,c);upload(directions,a.directions,c);
    upload(hits,a.hits,c);upload(summaries,a.summaries,c);upload(offsets,a.offsets,c);
    upload(lookup,hm,c);upload(count,hc,c);upload(records,hr,c);upload(result,ho,c);upload(focal,hf,c);
    FoldedPrepareParams prep{};prep.vertices=v.data();prep.patches=patches.data();
    prep.vertex_count=static_cast<unsigned>(a.vertices.size());prep.patch_count=static_cast<unsigned>(a.patches.size());
    prep.record_indices=lookup.data()+1;prep.records=records.data()+1;prep.written_count=count.data()+1;
    prep.capacity=capacity;prep.incident_direction=a.incident;prep.with_focal=1;prep.focal=a.focal;
    launch(module,"prepare_folded_patches",prep,prep.patch_count,c);
    FoldedOpticsParams q{};q.wave.optical=a.view();q.wave.optical.vertices=v.data();q.wave.optical.patches=patches.data();
    q.wave.optical.directions=directions.data();q.wave.optical.hits=hits.data();q.wave.optical.offsets=offsets.data();
    q.wave.optical.summaries=summaries.data();q.wave.optical.results=result.data()+1;
    q.wave.focal=a.focal;q.wave.focal_results=focal.data()+1;q.with_focal=1;
    q.folded={lookup.data()+1,records.data()+1,prep.patch_count,capacity};
    launch(module,"evaluate_folded_optics",q,n,c,128);
    RAINBOW_CUDA_CHECK(cuStreamSynchronize(c.stream()));
    lookup.download(hm);count.download(hc);records.download(hr);result.download(ho);focal.download(hf);
    require(hc[1]==capacity&&hc.front()==0xabcdef01u&&hc.back()==0xabcdef01u,"GPU folded counter guard failed.");
    require(hm.front()==0xabcdef01u&&hm.back()==0xabcdef01u&&hr.front().flags==0xabcdef01u&&hr.back().flags==0xabcdef01u,"GPU folded records guard failed.");
    require(ho.front().reserved==0xabcdef01u&&ho.back().reserved==0xabcdef01u&&hf.front().flags==0xabcdef01u&&hf.back().flags==0xabcdef01u,"GPU output guard failed.");
    for(const auto& r:host.records)
    {
        const auto idx=hm[r.compact_index+1];require(idx<capacity,"Missing GPU folded record.");
        const auto& actual=hr[idx+1];check_record(actual);
        require(actual.patch_id==r.patch_id,"GPU record index mismatch.");
        for(unsigned k=0;k<r.branch_count;++k)
            require(near(actual.branches[k].area_fraction,r.branches[k].area_fraction,1e-8)
                &&near(actual.branches[k].solid_angle_sr,r.branches[k].solid_angle_sr,1e-8,1e-16),"GPU branch metric mismatch.");
    }
    for(unsigned i=0;i<n;++i)
    {
        FocalOpticalResult reference_f{};
        const auto reference=PatchOpticalEvaluator::evaluate_direction(a.view(),i,&a.focal,&reference_f,host.view());
        check_result_pair(ho[i+1],hf[i+1],reference,reference_f);
        require(ho[i+1].known_hits_complete()&&hf[i+1].valid(),"GPU witness is still pending.");
    }
    std::cout<<"GPU snapshot replay: all 36 directions evaluated, including 12 recovered contributions.\n";
}
void production_api(const CudaContext& c,const std::filesystem::path& optical,
                    const std::filesystem::path& build,const std::filesystem::path& query_path)
{
    RaindropTraceConfig config{};config.grid_width=config.grid_height=2;config.grid_half_extent=1;
    config.incident_direction={0,0,1};config.incident_basis_x={1,0,0};config.radius_mm=.4f;
    Vec3 corners[4]={{0,0,1},Vec3{.2f,0,1}.normalized(),Vec3{.2f,0,1}.normalized(),Vec3{.4f,.1f,1}.normalized()};
    std::vector<OutgoingVertex> host(16);
    for(unsigned f=0;f<4;++f)for(unsigned j=0;j<4;++j)
    {
        auto& v=host[4*f+j];v.status=VertexStatus::Valid;v.direction_drop=corners[j];
        v.basis_x=TransverseFrame::from_direction(v.direction_drop).e0;v.field={{1,0},{0,0}};v.optical_cycles={1000,0};
    }
    OptixContext optix(c);DeviceBuffer<OutgoingVertex> vertices(c);upload(vertices,host,c);
    PatchAccel accel(c,optix);accel.load_module(build);accel.build(vertices,config);
    PatchQuery query(c,optix);query.create_pipeline(query_path);query.query_grid(accel,config,90,180);
    PatchOptics optics(c);optics.load_module(optical);optics.evaluate_wave(accel,query,config);
    require(optics.folded_statistics().prepared==4,"Production API did not prepare four folded primitives.");
    require(optics.folded_statistics().evaluated_hits>0,"Production API evaluated no branch contributions.");
    require(optics.statistics().error_directions==0,"Production folded API reported numeric errors.");
    auto first=std::vector<PatchOpticalResult>(optics.host_results().begin(),optics.host_results().end());
    optics.evaluate_wave(accel,query,config);
    require(optics.folded_statistics().prepared==4,"Repeated preparation lost candidates.");
    for(std::size_t i=0;i<first.size();++i)
        require(first[i].flags==optics.host_results()[i].flags,"Repeat changed availability.");
    optics.enable_folded_patches(false);optics.evaluate_wave(accel,query,config);
    require(optics.folded_statistics().evaluated_hits==0&&optics.statistics().pending_directions>0,"Legacy mode does not retain original pending behavior.");
    optics.enable_folded_patches(true);optics.evaluate(accel,query,config);
    require(!optics.has_wave_result()&&optics.folded_statistics().evaluated_hits>0,"Path-only production API does not use the branch model.");
    optics.close();query.close();accel.close();vertices.close();
    std::cout<<"Production API: original GAS -> query -> branch preparation -> optics; repeat/legacy/path-only passed.\n";
}
}
#if defined(_WIN32)
int wmain(int argc,wchar_t** argv)
#else
int main(int argc,char** argv)
#endif
{
    try
    {
        if(argc!=4)throw std::invalid_argument("Usage: folded_cuda <patch_optics.fatbin> <patch_build.fatbin> <patch_query.optixir>");
        CudaContext c{0};witness(c,argv[1]);production_api(c,argv[1],argv[2],argv[3]);
        std::cout<<folded_tests::checks<<" folded-patch checks passed (GPU).\n";return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
