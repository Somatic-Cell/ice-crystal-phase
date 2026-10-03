#include "patch_optics_test_cases.hpp"
#include "patch_query_validation.hpp"
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
using namespace rainbow;

// Test-only ownership for uploading the independent arithmetic fixtures. The
// production kernel is loaded unchanged; no duplicate arithmetic implementation.
class FixtureBuffers final
{
public:
    explicit FixtureBuffers(const CudaContext& cuda) noexcept
        : cuda_(cuda), module_(cuda), vertices_(cuda), patches_(cuda), directions_(cuda),
          offsets_(cuda), summaries_(cuda), hits_(cuda), results_(cuda) {}
    ~FixtureBuffers() noexcept
    {
        if(pending_)
        {
            static_cast<void>(detail::report_cuda_cleanup_result(cuCtxSetCurrent(cuda_.handle()),"test context"));
            static_cast<void>(detail::report_cuda_cleanup_result(cuStreamSynchronize(cuda_.stream()),"test sync"));
        }
    }
    FixtureBuffers(const FixtureBuffers&) = delete;
    FixtureBuffers& operator=(const FixtureBuffers&) = delete;

    void run(const std::filesystem::path& path, const tests::OpticsFixture& input)
    {
        cuda_.make_current();module_.load_fatbin(path);
        vertices_.allocate(input.vertices.size());patches_.allocate(input.patches.size());
        directions_.allocate(input.directions.size());offsets_.allocate(input.offsets.size());
        summaries_.allocate(input.summaries.size());hits_.allocate(input.hits.size());
        // Guard records verify the last partially occupied CUDA block.
        host_output_.assign(input.directions.size()+2u,PatchOpticalResult{});
        host_output_.front().reserved=0x12345678u;host_output_.back().reserved=0xabcdef12u;
        results_.allocate(host_output_.size());
        pending_=true;
        vertices_.upload_async(input.vertices,cuda_.stream());patches_.upload_async(input.patches,cuda_.stream());
        directions_.upload_async(input.directions,cuda_.stream());offsets_.upload_async(input.offsets,cuda_.stream());
        summaries_.upload_async(input.summaries,cuda_.stream());hits_.upload_async(input.hits,cuda_.stream());
        results_.upload_async(host_output_,cuda_.stream());
        auto params=input.view();
        params.vertices=vertices_.data();params.patches=patches_.data();params.directions=directions_.data();
        params.offsets=offsets_.data();params.summaries=summaries_.data();params.hits=hits_.data();
        params.results=results_.data()+1;
        void* args[]={&params};
        constexpr unsigned threads=128;
        const unsigned blocks=params.direction_count/threads+(params.direction_count%threads!=0u ? 1u:0u);
        RAINBOW_CUDA_CHECK(cuLaunchKernel(module_.find_function("evaluate_patch_optics"),
            blocks,1,1,threads,1,1,0,cuda_.stream(),args,nullptr));
        RAINBOW_CUDA_CHECK(cuStreamSynchronize(cuda_.stream()));pending_=false;
        results_.download(host_output_);
        tests::require_optics(host_output_.front().reserved==0x12345678u
            &&host_output_.back().reserved==0xabcdef12u,"CUDA output overrun.");
        const std::span<const PatchOpticalResult> actual{host_output_.data()+1,input.directions.size()};
        tests::verify_optics_fixture(input,actual);
        for(std::uint32_t i=0;i<params.direction_count;++i)
            tests::compare_optics_results(actual[i],PatchOpticalEvaluator::evaluate_direction(input.view(),i));
        results_.close();hits_.close();summaries_.close();offsets_.close();directions_.close();
        patches_.close();vertices_.close();module_.close();
        std::cout<<"CUDA optical arithmetic: "<<actual.size()<<" cases passed.\n";
    }
private:
    const CudaContext& cuda_;
    CudaModule module_;
    DeviceBuffer<OutgoingVertex> vertices_;
    DeviceBuffer<OutgoingPatch> patches_;
    DeviceBuffer<Vec3> directions_;
    DeviceBuffer<std::uint64_t> offsets_;
    DeviceBuffer<PatchQuerySummary> summaries_;
    DeviceBuffer<PatchQueryHit> hits_;
    DeviceBuffer<PatchOpticalResult> results_;
    std::vector<PatchOpticalResult> host_output_;
    bool pending_=false;
};

int run(const std::filesystem::path& optics_module, const std::filesystem::path& build_module,
        const std::filesystem::path& query_module)
{
    CudaContext cuda{0};
    const auto fixture=tests::make_optics_fixture(); // must outlive asynchronous uploads
    { FixtureBuffers buffers(cuda);buffers.run(optics_module,fixture); }

    // Actual CUDA patch construction -> OptiX query -> PatchOptics production API.
    // Reuse the existing curved/folded geometry fixture instead of a new solver.
    OptixContext optix(cuda);
    auto host_vertices=tests::query_test_vertices();
    for(std::size_t i=0;i<host_vertices.size();++i)
        host_vertices[i].optical_cycles={1000+static_cast<std::int32_t>(i/4),float(i%4)/8.0f};
    DeviceBuffer<OutgoingVertex> vertices(cuda);
    vertices.allocate(host_vertices.size());
    vertices.upload_async(host_vertices,cuda.stream());
    RAINBOW_CUDA_CHECK(cuStreamSynchronize(cuda.stream()));
    RaindropTraceConfig config{};config.grid_width=2;config.grid_height=2;
    PatchAccel accel(cuda,optix);accel.load_module(build_module);accel.build(vertices,config);
    PatchQuery query(cuda,optix);query.create_pipeline(query_module);
    const auto g=tests::front_patch();
    const std::vector<Vec3> directions={{0,0,1},{0,0,-1},{1,0,0},g.corners[0],
        Vec3{.11f,.05f,1}.normalized(),Vec3{-.1f,.03f,1}.normalized()};
    query.query(accel,directions);
    tests::require_optics(query.statistics().error_directions==0u,"Query fixture failed before optics.");
    PatchOptics optics(cuda);optics.load_module(optics_module);optics.evaluate(accel,query,config);
    tests::require_optics(optics.statistics().error_directions==0u,"Integration optical numerical error.");
    tests::require_optics(optics.statistics().pending_directions!=0u,"Folded/boundary hits were silently treated as complete.");
    std::vector<OutgoingPatch> host_patches(accel.patches().element_count());
    accel.patches().download(host_patches);
    const auto host_hits=query.download_hits();
    PatchOpticsParams reference{};
    reference.vertices=host_vertices.data();reference.patches=host_patches.data();
    reference.directions=directions.data();reference.offsets=query.host_offsets().data();
    reference.summaries=query.host_summaries().data();reference.hits=host_hits.data();
    reference.vertex_count=static_cast<std::uint32_t>(host_vertices.size());
    reference.patch_count=static_cast<std::uint32_t>(host_patches.size());
    reference.direction_count=static_cast<std::uint32_t>(directions.size());
    reference.hit_storage_count=host_hits.size();reference.incident_direction=config.incident_direction;
    reference.incident_basis_x=config.incident_basis_x;
    for(std::uint32_t i=0;i<reference.direction_count;++i)
        tests::compare_optics_results(optics.host_results()[i],PatchOpticalEvaluator::evaluate_direction(reference,i));
    const std::vector<PatchOpticalResult> first(optics.host_results().begin(),optics.host_results().end());
    optics.evaluate(accel,query,config);
    for(std::size_t i=0;i<first.size();++i) tests::compare_optics_results(optics.host_results()[i],first[i]);
    // Empty GAS, then empty direction list. No fabricated optical failures/NaNs.
    for(auto& v:host_vertices)v.status=VertexStatus::Miss;
    vertices.upload_async(host_vertices,cuda.stream());
    RAINBOW_CUDA_CHECK(cuStreamSynchronize(cuda.stream()));
    accel.build(vertices,config);query.query(accel,directions);optics.evaluate(accel,query,config);
    tests::require_optics(optics.statistics().evaluated_hits==0u
        &&optics.statistics().no_outgoing_source_cells==4u,"Empty source metadata lost.");
    for(const auto& r:optics.host_results())tests::require_optics(r.known_hits_complete()
        &&r.regular_partial_path_s==0.0,"Empty known-hit set not handled.");
    query.query(accel,std::span<const Vec3>{});optics.evaluate(accel,query,config);
    tests::require_optics(optics.host_results().empty()&&optics.statistics().directions==0u,"Zero directions.");
    optics.close();query.close();accel.close();vertices.close();
    std::cout<<"Patch optics integration: CUDA build, OptiX hits, CPU/GPU comparison, repeat, empty passed.\n";
    return EXIT_SUCCESS;
}
} // namespace

#if defined(_WIN32)
int wmain(const int argc,wchar_t* argv[])
#else
int main(const int argc,char* argv[])
#endif
{
    try
    {
        if(argc!=4)throw std::invalid_argument(
            "Usage: optics_test <patch_optics.fatbin> <patch_build.fatbin> <patch_query.optixir>");
        return run(std::filesystem::path{argv[1]},std::filesystem::path{argv[2]},std::filesystem::path{argv[3]});
    }
    catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<'\n';return EXIT_FAILURE;}
}
