#include "unpolarized_test_support.hpp"
#include <rainbow/raindrop_tracer.hpp>
#include <rainbow/patch_optics.hpp>
#include <rainbow/patch_query.hpp>
#include <rainbow/rainbow_diffraction.hpp>
#include <rainbow/patch_failure_report.hpp>
#include <chrono>
#include <fstream>
#include <iterator>
#include <string>
#include <iomanip>
#include <span>
#include <exception>

using namespace unpolarized_tests;
namespace
{
std::string read_closed(const std::filesystem::path& p)
{
    std::ifstream file(p,std::ios::binary);
    if(!file)throw std::runtime_error("Cannot read test output.");
    return {std::istreambuf_iterator<char>(file),std::istreambuf_iterator<char>()};
}
struct TestDirectory
{
    std::filesystem::path path;
    TestDirectory()
    {
        path=std::filesystem::temp_directory_path()/
            ("rainbow_unpolarized_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        if(!std::filesystem::create_directory(path))throw std::runtime_error("Test output directory already exists.");
    }
    ~TestDirectory(){std::error_code error;std::filesystem::remove_all(path,error);}
};
// Failure-only diagnostics. These functions do not change the pass/fail
// criterion, production buffers, quadrature settings, or optical algorithms.
// Three replays use the existing production CUDA kernel and identical stored
// vertices/hits/CPU-prepared folded records; no ray tracing is repeated here.
struct OpticalReplay
{
    PatchOpticalResult optical{};
    FocalOpticalResult focal{};
};
struct GpuReplaySet
{
    OpticalReplay paired{}, input0{}, input1{};
};

template<class T>
void upload_snapshot(DeviceBuffer<T>& destination, std::span<const T> source)
{
    destination.allocate(source.size());
    if(!source.empty())
        RAINBOW_CUDA_CHECK(cuMemcpyHtoD(
            destination.address(), source.data(), destination.byte_size()));
}

GpuReplaySet replay_on_gpu(CudaContext& context,
                          const std::filesystem::path& module_path,
                          const Traced& paired, const Traced& input1,
                          const Geometry& geometry, const Query& query)
{
    context.make_current();
    CudaModule module(context);
    module.load_fatbin(module_path);
    static_cast<void>(module.find_function("rainbow_two_input_optics_abi_v2"));
    const auto kernel = module.find_function("evaluate_folded_optics");

    DeviceBuffer<OutgoingVertex> vertices0(context), vertices1(context);
    DeviceBuffer<Field32> second(context);
    DeviceBuffer<OutgoingPatch> patches(context);
    DeviceBuffer<std::uint32_t> lookup(context);
    DeviceBuffer<FoldedPatchRecord> folded(context);
    DeviceBuffer<Vec3> directions(context);
    DeviceBuffer<std::uint64_t> offsets(context);
    DeviceBuffer<PatchQuerySummary> summaries(context);
    DeviceBuffer<PatchQueryHit> hits(context);
    DeviceBuffer<PatchOpticalResult> result(context);
    DeviceBuffer<FocalOpticalResult> focal_result(context);
    // Deliberately use synchronous uploads. No input may outlive an unfinished
    // transfer on an allocation or kernel error path.
    upload_snapshot(vertices0, std::span<const OutgoingVertex>(paired.vertices));
    upload_snapshot(vertices1, std::span<const OutgoingVertex>(input1.vertices));
    upload_snapshot(second, std::span<const Field32>(paired.second));
    upload_snapshot(patches, std::span<const OutgoingPatch>(geometry.patches));
    upload_snapshot(lookup, std::span<const std::uint32_t>(geometry.lookup));
    upload_snapshot(folded, std::span<const FoldedPatchRecord>(geometry.folded));
    upload_snapshot(directions, std::span<const Vec3>(&query.direction, 1));
    upload_snapshot(offsets, std::span<const std::uint64_t>(query.offsets, 2));
    upload_snapshot(summaries, std::span<const PatchQuerySummary>(&query.summary, 1));
    upload_snapshot(hits, std::span<const PatchQueryHit>(query.hits));
    result.allocate(1);
    focal_result.allocate(1);

    const auto launch = [&](bool dual, bool use_input1)
    {
        FoldedOpticsParams parameters{};
        auto& optical = parameters.wave.optical;
        optical = query.params(paired, geometry, dual);
        optical.vertices = use_input1 ? vertices1.data() : vertices0.data();
        optical.patches = patches.data();
        optical.directions = directions.data();
        optical.offsets = offsets.data();
        optical.summaries = summaries.data();
        optical.hits = hits.data();
        optical.results = result.data();
        optical.second_input_fields = dual ? second.data() : nullptr;
        optical.second_input_count = dual ? optical.vertex_count : 0u;
        parameters.wave.focal = geometry.focal;
        parameters.wave.focal_results = focal_result.data();
        parameters.folded = {lookup.data(), folded.data(),
            static_cast<std::uint32_t>(geometry.lookup.size()),
            static_cast<std::uint32_t>(geometry.folded.size())};
        parameters.with_focal = 1u;
        void* arguments[] = {&parameters};
        try
        {
            RAINBOW_CUDA_CHECK(cuLaunchKernel(kernel, 1, 1, 1, 128, 1, 1, 0,
                context.stream(), arguments, nullptr));
            RAINBOW_CUDA_CHECK(cuStreamSynchronize(context.stream()));
        }
        catch(...)
        {
            // Keep buffers and the module alive until the completion attempt.
            const auto status = cuStreamSynchronize(context.stream());
            if(status != CUDA_SUCCESS)
                std::cerr << "[u2diag] replay cleanup synchronization status="
                          << static_cast<int>(status) << '\n';
            throw;
        }
        OpticalReplay output;
        result.download(std::span<PatchOpticalResult>(&output.optical, 1));
        focal_result.download(std::span<FocalOpticalResult>(&output.focal, 1));
        return output;
    };
    GpuReplaySet output;
    output.paired = launch(true, false);
    output.input0 = launch(false, false);
    output.input1 = launch(false, true);
    return output;
}

void compare_scalar(const std::string& name, double a, double b, double tolerance)
{
    const double error = std::abs(a - b);
    const double limit = 1e-12 + tolerance * (std::max)(std::abs(a), std::abs(b));
    std::cerr << "[u2diag] " << name << " actual=" << a << " reference=" << b
              << " abs_error=" << error << " allowed=" << limit
              << " verdict=" << (near(a,b,tolerance) ? "PASS" : "FAIL") << '\n';
}

void compare_field(const std::string& name, const OpticalField64& a,
                   const OpticalField64& b, double tolerance)
{
    compare_scalar(name+".s_real",a.s_real,b.s_real,tolerance);
    compare_scalar(name+".s_imag",a.s_imag,b.s_imag,tolerance);
    compare_scalar(name+".p_real",a.p_real,b.p_real,tolerance);
    compare_scalar(name+".p_imag",a.p_imag,b.p_imag,tolerance);
    const double error = std::hypot(std::hypot(a.s_real-b.s_real,a.s_imag-b.s_imag),
                                    std::hypot(a.p_real-b.p_real,a.p_imag-b.p_imag));
    const double scale = (std::max)(
        std::hypot(std::hypot(a.s_real,a.s_imag),std::hypot(a.p_real,a.p_imag)),
        std::hypot(std::hypot(b.s_real,b.s_imag),std::hypot(b.p_real,b.p_imag)));
    // This is an additional diagnostic only. It is NOT a replacement threshold.
    std::cerr << "[u2diag] " << name << " field_norm=" << scale
              << " error_norm=" << error << " norm_relative_error="
              << (scale>0 ? error/scale : error) << '\n';
}

void describe_average(const char* name, const OpticalReplay& paired,
                      const OpticalReplay& x, const OpticalReplay& y, double tolerance)
{
    try
    {
        check_average(paired.optical, paired.focal, x.optical, x.focal,
                      y.optical, y.focal, tolerance);
        std::cerr << "[u2diag] " << name << " PASS\n";
    }
    catch(const std::exception& error)
    {
        // Only the diagnostic sub-check is caught. The original failure is
        // rethrown by the caller, so this cannot turn a failed CTest into a pass.
        std::cerr << "[u2diag] " << name << " FAIL: " << error.what() << '\n';
    }
}

void diagnose_optical_mismatch(CudaContext& context,
    const std::filesystem::path& module_path, const RaindropSettings& settings,
    std::uint32_t direction_id, std::uint32_t theta_count, std::uint32_t phi_count,
    const Traced& paired, const Traced& y, const Geometry& geometry, const Query& query,
    const OpticalReplay& actual, const OpticalReplay& host_x, const OpticalReplay& host_y)
{
    constexpr double tolerance = 4e-5; // Exactly the original test threshold.
    const auto old_precision = std::cerr.precision();
    std::cerr << std::setprecision(17)
              << "\n[u2diag] BEGIN: original assertion remains a failure\n"
              << "[u2diag] radius_mm=" << settings.radius_mm
              << " force_sphere=" << settings.force_sphere
              << " grid=" << settings.grid_width << 'x' << settings.grid_height
              << " direction_id=" << direction_id
              << " theta_deg=" << 180.0*(double(direction_id/phi_count)+0.5)/theta_count
              << " phi_deg=" << -180.0+360.0*(double(direction_id%phi_count)+0.5)/phi_count << '\n'
              << "[u2diag] query_flags=" << query.summary.flags
              << " optical_flags=" << actual.optical.flags
              << " focal_flags=" << actual.focal.flags
              << " hit_count=" << actual.optical.hit_count
              << " evaluated_hits=" << actual.optical.evaluated_hits
              << " rejected_hits=" << actual.optical.rejected_hits << '\n';
    compare_field("original_GPU_vs_CPU.input0.path", actual.optical.regular_partial_path_field,
                  host_x.optical.regular_partial_path_field,tolerance);
    compare_field("original_GPU_vs_CPU.input1.path", actual.optical.regular_partial_path_field_second,
                  host_y.optical.regular_partial_path_field,tolerance);
    compare_field("original_GPU_vs_CPU.input1.focal",actual.focal.field_second,host_y.focal.field,tolerance);
    compare_scalar("unpolarized.incoherent_s",actual.optical.regular_partial_incoherent_s,
        0.5*(host_x.optical.regular_partial_incoherent_s+host_y.optical.regular_partial_incoherent_s),tolerance);
    compare_scalar("unpolarized.incoherent_p",actual.optical.regular_partial_incoherent_p,
        0.5*(host_x.optical.regular_partial_incoherent_p+host_y.optical.regular_partial_incoherent_p),tolerance);
    compare_scalar("unpolarized.path_s",actual.optical.regular_partial_path_s,
        0.5*(host_x.optical.regular_partial_path_s+host_y.optical.regular_partial_path_s),tolerance);
    compare_scalar("unpolarized.path_p",actual.optical.regular_partial_path_p,
        0.5*(host_x.optical.regular_partial_path_p+host_y.optical.regular_partial_path_p),tolerance);
    compare_scalar("unpolarized.focal_s",actual.focal.intensity_s,
        0.5*(host_x.focal.intensity_s+host_y.focal.intensity_s),tolerance);
    compare_scalar("unpolarized.focal_p",actual.focal.intensity_p,
        0.5*(host_x.focal.intensity_p+host_y.focal.intensity_p),tolerance);
    for(std::uint32_t i=0; i<query.summary.hit_count && i<query.hits.size(); ++i)
    {
        const auto& h = query.hits[i];
        std::cerr << "[u2diag] hit=" << i << " patch_id=" << h.patch_id
                  << " compact_index=" << h.compact_index << " root_index=" << h.root_index
                  << " flags=" << h.flags << " u=" << h.u << " v=" << h.v
                  << " t=" << h.t << " residual=" << h.residual;
        if(h.compact_index < geometry.patches.size())
        {
            const auto& patch = geometry.patches[h.compact_index];
            std::cerr << " stored_patch_status=" << static_cast<std::uint32_t>(patch.status);
        }
        std::cerr << '\n';
    }
    try
    {
        OpticalReplay host_pair;
        host_pair.optical = PatchOpticalEvaluator::evaluate_direction(
            query.params(paired,geometry,true),0,&geometry.focal,&host_pair.focal,geometry.view());
        describe_average("CPU_paired_vs_CPU_singles",host_pair,host_x,host_y,tolerance);
        const auto gpu = replay_on_gpu(context,module_path,paired,y,geometry,query);
        std::cerr << "[u2diag] GPU replay: same stored inputs and CPU-prepared folded records; "
                     "existing evaluate_folded_optics kernel; no new ray trace\n";
        describe_average("GPU_paired_vs_GPU_singles",gpu.paired,gpu.input0,gpu.input1,tolerance);
        compare_field("original_GPU_vs_replay_GPU.input0.path",actual.optical.regular_partial_path_field,
                      gpu.paired.optical.regular_partial_path_field,tolerance);
        compare_field("original_GPU_vs_replay_GPU.input1.path",actual.optical.regular_partial_path_field_second,
                      gpu.paired.optical.regular_partial_path_field_second,tolerance);
        compare_field("replay_GPU_vs_CPU.input1.path",gpu.input1.optical.regular_partial_path_field,
                      host_y.optical.regular_partial_path_field,tolerance);
        compare_field("replay_GPU_pair_vs_GPU_single.input1.path",gpu.paired.optical.regular_partial_path_field_second,
                      gpu.input1.optical.regular_partial_path_field,tolerance);
    }
    catch(const std::exception& error)
    {
        std::cerr << "[u2diag] diagnostic replay itself failed: " << error.what() << '\n';
    }
    std::cerr << "[u2diag] END: original test failure is rethrown\n\n";
    std::cerr.precision(old_precision);
}

void run_case(CudaContext& c,RaindropTracer& tracer,const RaindropSettings& s,
              const std::filesystem::path& patch_module,const std::filesystem::path& query_module,
              const std::filesystem::path& optics_module)
{
    // Test-only separate inputs establish a reference. Production does not retrace.
    tracer.trace_unpolarized(s);
    Traced paired;paired.config=tracer.config();paired.vertices=tracer.download_vertices();paired.second=tracer.download_second_input_fields();
    auto sx=s,sy=s;sx.incident_field={{1,0},{0,0}};sy.incident_field={{0,0},{1,0}};
    tracer.trace(sx);const auto vx=tracer.download_vertices();
    require(!tracer.is_unpolarized()&&tracer.second_input_fields().is_empty(),"Legacy retrace did not clear paired state.");
    tracer.trace(sy);const auto vy=tracer.download_vertices();
    for(std::size_t i=0;i<vx.size();++i)
    {
        require(same_geometry(paired.vertices[i],vx[i])&&same_geometry(paired.vertices[i],vy[i]),"GPU paired trace changed geometry.");
        if(vx[i].status==VertexStatus::Valid)
        {
            require(near(paired.vertices[i].field,vx[i].field,1e-5),"GPU input-0 field differs from single-input run.");
            require(near(paired.second[i],vy[i].field,1e-5),"GPU input-1 field differs from single-input run.");
        }
    }
    tracer.trace_unpolarized(s); // new source, before ANY GAS borrows it
    paired.vertices=tracer.download_vertices();paired.second=tracer.download_second_input_fields();
    PatchAccel accel(c,tracer.optix_context());accel.load_module(patch_module);accel.build(tracer);
    require(accel.is_unpolarized(),"Paired state was lost during GAS build.");
    PatchQuery query(c,tracer.optix_context());query.create_pipeline(query_module);
    constexpr unsigned nt=90,np=180;
    query.query_grid(accel,tracer.config(),nt,np);
    PatchOptics optics(c);optics.load_module(optics_module);optics.evaluate_wave(accel,query,tracer.config());
    require(optics.is_unpolarized()&&optics.has_wave_result(),"Optics did not publish unpolarized data.");
    Geometry geometry;
    geometry.patches.resize(accel.patches().element_count());accel.patches().download(geometry.patches);
    geometry.lookup.assign(geometry.patches.size(),0xffffffffu);
    geometry.focal.grid_width=s.grid_width;geometry.focal.grid_height=s.grid_height;
    geometry.focal.grid_half_extent=paired.config.grid_half_extent;
    for(unsigned i=0;i<geometry.patches.size();++i)
    {
        if(geometry.patches[i].status!=PatchCellStatus::NeedsRefinement)continue;
        geometry.lookup[i]=static_cast<unsigned>(geometry.folded.size());
        geometry.folded.push_back(FoldedPatchBuilder::prepare(paired.vertices.data(),static_cast<unsigned>(paired.vertices.size()),
            geometry.patches[i],i,paired.config.incident_direction,&geometry.focal));
    }
    Traced x=paired,y=paired;x.second.clear();y.second.clear();
    for(std::size_t i=0;i<y.vertices.size();++i)y.vertices[i].field=paired.second[i];
    const auto hits=query.download_hits();
    const auto& offsets=query.host_offsets();const auto& summaries=query.host_summaries();const auto& directions=query.host_directions();
    std::vector<FocalOpticalResult> expected_focal(directions.size());
    unsigned complete=0;
    for(unsigned i=0;i<directions.size();++i)
    {
        Query q;q.direction=directions[i];q.summary=summaries[i];
        q.hits.assign(hits.begin()+static_cast<std::ptrdiff_t>(offsets[i]),hits.begin()+static_cast<std::ptrdiff_t>(offsets[i+1]));
        q.offsets[1]=q.hits.size();
        FocalOpticalResult fx{},fy{};
        const auto px=PatchOpticalEvaluator::evaluate_direction(q.params(x,geometry,false),0,&geometry.focal,&fx,geometry.view());
        const auto py=PatchOpticalEvaluator::evaluate_direction(q.params(y,geometry,false),0,&geometry.focal,&fy,geometry.view());
        try
        {
            check_average(optics.host_results()[i],optics.host_focal_results()[i],px,fx,py,fy,4e-5);
        }
        catch(const std::exception&)
        {
            const auto original_failure = std::current_exception();
            try
            {
                diagnose_optical_mismatch(c,optics_module,s,i,nt,np,paired,y,geometry,q,
                    {optics.host_results()[i],optics.host_focal_results()[i]},
                    {px,fx},{py,fy});
            }
            catch(const std::exception& diagnostic_error)
            {
                std::cerr << "[u2diag] diagnostic exception: " << diagnostic_error.what() << '\n';
            }
            catch(...)
            {
                std::cerr << "[u2diag] unexpected diagnostic exception\n";
            }
            std::rethrow_exception(original_failure); // Preserve the original failed assertion.
        }
        auto& f=expected_focal[i];f=fx;f.flags|=fy.flags;
        f.intensity_s=.5*(fx.intensity_s+fy.intensity_s);f.intensity_p=.5*(fx.intensity_p+fy.intensity_p);
        for(unsigned k=0;k<4;++k)f.family_incoherent[k]=.5*(fx.family_incoherent[k]+fy.family_incoherent[k]);
        complete+=f.valid()&&px.known_hits_complete();
    }
    require(complete>directions.size()/2,"Too few complete GPU optical comparisons.");
    std::vector<RainbowTransition> transitions(directions.size());
    DiffractionParams p{};p.focal=expected_focal.data();p.theta_count=nt;p.phi_count=np;p.transitions=transitions.data();
    double sigma=0;require(RainbowDiffraction::table_sigma_degrees(s.radius_mm,sigma),"Missing bandwidth.");
    p.config.primary_sigma_rad=sigma*RainbowDiffraction::pi/180;
    for(unsigned i=0;i<directions.size();++i)transitions[i]=RainbowDiffraction::detect(p,i);
    for(unsigned i=0;i<directions.size();++i)
    {
        const auto expected=RainbowDiffraction::filter(p,i);const auto actual=optics.host_diffraction_results()[i];
        require(actual.flags==expected.flags,"Diffraction was not detected after input averaging.");
        if(expected.valid())require(near(expected.intensity_s,actual.intensity_s,5e-5)&&near(expected.intensity_p,actual.intensity_p,5e-5),"Diffraction differs from filter(mean of two focal inputs).");
    }
    TestDirectory dir;
    tracer.write_csv(dir.path/"vertices.csv");optics.write_csv(dir.path/"optics.csv");optics.write_wave_csv(dir.path/"wave.csv");
    require(read_closed(dir.path/"vertices.csv").find("rainbow_outgoing_vertices_v2")!=std::string::npos,"Vertex CSV version not updated.");
    const auto optical_text=read_closed(dir.path/"optics.csv"),wave_text=read_closed(dir.path/"wave.csv");
    require(optical_text.find("# incident_polarization=unpolarized")!=std::string::npos&&optical_text.find("path_J11_imag")!=std::string::npos,"Optical CSV lacks paired metadata/column.");
    require(wave_text.find("rainbow_wave_optics_v2")!=std::string::npos&&wave_text.find("focal_J11_imag")!=std::string::npos,"Wave CSV lacks paired metadata/column.");
    auto report=PatchFailureReport::capture(c,accel,query,optics,tracer.config(),{});report.write_json(dir.path/"report.json");
    require(report.snapshot().input_polarization==IncidentPolarization::Unpolarized,"Diagnostic capture lost input semantics.");
    require(read_closed(dir.path/"report.json").find("rainbow_patch_failure_witness_v2")!=std::string::npos,"Diagnostic JSON version not updated.");
    optics.close();query.close();accel.close();
    std::cout<<"GPU: "<<s.radius_mm<<" mm, grid "<<s.grid_width<<", "<<complete<<" complete directions; paired trace, mean-before-diffraction, CSV/JSON passed.\n";
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
        if(argc!=5)throw std::invalid_argument("Usage: unpolarized_cuda <raindrop.optixir> <patch_build.fatbin> <patch_query.optixir> <patch_optics.fatbin>");
        CudaContext c{0};RaindropTracer tracer(c);tracer.create_pipeline(argv[1]);
        run_case(c,tracer,settings(.4f,17,true),argv[2],argv[3],argv[4]);
        run_case(c,tracer,settings(1.0f,129),argv[2],argv[3],argv[4]);
        tracer.close();std::cout<<checks<<" unpolarized checks passed (GPU).\n";return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
