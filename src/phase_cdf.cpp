#include <rainbow/phase_cdf.hpp>
#include <rainbow/cuda_host_upload.hpp>
#include <rainbow/phase_cdf_math.hpp>
#include <rainbow/npy_writer.hpp>
#include <rainbow/patch_optics.hpp>
#include <rainbow/raindrop_tracer.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace rainbow
{
namespace
{
void validate_policy(const PhaseCdfPolicy& p)
{
    if(!(std::isfinite(p.maximum_l1_error) && p.maximum_l1_error>=0 && p.maximum_l1_error<=2 &&
         std::isfinite(p.maximum_lost_mass) && p.maximum_lost_mass>=0 && p.maximum_lost_mass<=1))
        throw std::invalid_argument("Invalid CDF rounding-error policy.");
}
std::string describe(const PhaseCdfReduction& s)
{
    std::ostringstream out;
    out << "max_intensity=" << s.maximum << ", invalid_values=" << s.invalid_values << ", incomplete=" << s.incomplete
        << ", underresolved=" << s.underresolved << ", malformed=" << s.malformed
        << ", first_problem=" << s.first_problem << ", L1=" << s.l1_error
        << ", lost_mass=" << s.lost_mass << ", weight_underflow=" << s.weight_underflow;
    return out.str();
}
std::string describe_moment(const PhaseMomentSum& value, const PhaseStorageParams& p)
{
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(17)
        << " mass=" << value.mass << " axial=" << value.axial
        << " weight_underflow=" << value.weight_underflow
        << " grid=" << p.nt << 'x' << p.np << " count=" << p.input.count
        << " divisor=" << p.divisor;
    return out.str();
}
// Failure-only readback of coordinate metadata, not of the optical table.
// A diagnostic failure must never hide the original failed moment gate.
std::string inspect_moment_edges(const CudaContext& context,
    const PhaseStorageParams& p, std::span<const double> expected)
{
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(17);
    try
    {
        if(!p.fine_edges || expected.size()!=std::size_t(p.nt)+1)
            return " edge_inspection=unavailable";
        context.make_current();
        std::vector<double> actual(expected.size());
        RAINBOW_CUDA_CHECK(cuMemcpyDtoH(actual.data(),
            reinterpret_cast<CUdeviceptr>(p.fine_edges), actual.size()*sizeof(double)));
        std::size_t mismatches=0, invalid_intervals=0, first=actual.size();
        for(std::size_t i=0;i<actual.size();++i)
        {
            if(actual[i]!=expected[i])
            {
                ++mismatches;
                if(first==actual.size()) first=i;
            }
            if(i && !(std::isfinite(actual[i-1]) && std::isfinite(actual[i])
                      && actual[i]>actual[i-1])) ++invalid_intervals;
        }
        out << " edge_mismatches=" << mismatches
            << " invalid_edge_intervals=" << invalid_intervals;
        if(first<actual.size())
            out << " first_edge=" << first << " expected=" << expected[first]
                << " observed=" << actual[first];
    }
    catch(const std::exception& e)
    {
        out << " edge_inspection_failed=" << e.what();
    }
    return out.str();
}
const char* stage_name(PhaseDensityStage stage) noexcept
{
    switch(stage)
    {
    case PhaseDensityStage::Scalar: return "scalar";
    case PhaseDensityStage::Incoherent: return "incoherent";
    case PhaseDensityStage::Path: return "path";
    case PhaseDensityStage::Focal: return "focal";
    case PhaseDensityStage::Diffraction: return "diffraction";
    }
    return "unknown";
}

template<class T>
T read_failure_record(const CudaContext& context, const T* source,
                      std::uint64_t index, std::uint64_t count)
{
    if(!source || index>=count)
        throw std::out_of_range("Invalid failure-inspection source/index.");
    const auto base=reinterpret_cast<CUdeviceptr>(source);
    if(index>((std::numeric_limits<CUdeviceptr>::max)()-base)/sizeof(T))
        throw std::overflow_error("Failure-inspection address overflow.");
    T result{};
    context.make_current();
    RAINBOW_CUDA_CHECK(cuMemcpyDtoH(&result,base+index*sizeof(T),sizeof(T)));
    return result;
}

// These are snapshots of already computed output records, NOT a CPU replay.
// One record per nonempty stage is enough to locate its first hard failure.
// Counter totals describe all directions; the exemplars do not claim that all
// failures share the same flags or cause.
void append_failure_record(std::ostream& out, const CudaContext& context,
                           const PhaseDensityView& input, std::uint64_t index,
                           std::uint32_t nt, std::uint32_t np)
{
    const auto row=index/np, col=index%np;
    out << "[cdf-input] example_for=" << stage_name(input.stage)
        << " direction=" << index << " theta_index=" << row << " phi_index=" << col
        << " theta_deg=" << (double(row)+0.5)*(180.0/double(nt))
        << " phi_deg=" << -180.0+(double(col)+0.5)*(360.0/double(np)) << '\n';
    if(input.stage==PhaseDensityStage::Scalar)
    {
        out << "[cdf-input] scalar=" << read_failure_record(context,input.scalar,index,input.count) << '\n';
        return;
    }
    const auto o=read_failure_record(context,input.optical,index,input.count);
    out << "[cdf-input] optical flags=" << o.flags << " query_flags=" << o.query_flags
        << " known_hits_complete=" << o.known_hits_complete()
        << " hit_count=" << o.hit_count << " evaluated_hits=" << o.evaluated_hits
        << " rejected_hits=" << o.rejected_hits << " refinement_hits=" << o.refinement_hits
        << " boundary_hits=" << o.boundary_hits << " singular_hits=" << o.singular_hits
        << " folded_evaluated_hits=" << o.folded_evaluated_hits()
        << " first_problem_patch_id=" << o.first_problem_patch_id
        << " incoherent_s=" << o.regular_partial_incoherent_s
        << " incoherent_p=" << o.regular_partial_incoherent_p
        << " path_s=" << o.regular_partial_path_s
        << " path_p=" << o.regular_partial_path_p << '\n';
    if(input.stage==PhaseDensityStage::Focal || input.stage==PhaseDensityStage::Diffraction)
    {
        const auto f=read_failure_record(context,input.focal,index,input.count);
        out << "[cdf-input] focal flags=" << f.flags << " valid=" << f.valid()
            << " corrected_hits=" << f.corrected_hits
            << " first_problem_patch_id=" << f.first_problem_patch_id
            << " intensity_s=" << f.intensity_s << " intensity_p=" << f.intensity_p
            << " family_hits=" << f.family_hits[0] << ',' << f.family_hits[1]
            << ',' << f.family_hits[2] << ',' << f.family_hits[3] << '\n';
    }
    if(input.stage==PhaseDensityStage::Diffraction)
    {
        const auto d=read_failure_record(context,input.diffraction,index,input.count);
        out << "[cdf-input] diffraction flags=" << d.flags << " valid=" << d.valid()
            << " intensity_s=" << d.intensity_s << " intensity_p=" << d.intensity_p
            << " sigma_rad=" << d.sigma_rad << " blend=" << d.blend
            << " transition_kind=" << d.transition_kind << '\n';
    }
}

class PinnedDoubles final
{
public:
    explicit PinnedDoubles(std::size_t n) : count(n)
    { RAINBOW_CUDA_CHECK(cuMemHostAlloc(reinterpret_cast<void**>(&data), n*sizeof(double), 0)); }
    ~PinnedDoubles() noexcept
    { if(data) static_cast<void>(detail::report_cuda_cleanup_result(cuMemFreeHost(data), "cuMemFreeHost(CDF staging)")); }
    double* data=nullptr; std::size_t count;
};
void write_device_npy(const CudaContext& c, const DeviceBuffer<double>& buffer,
    const std::filesystem::path& path, std::span<const std::uint64_t> shape, PinnedDoubles& staging)
{
    NpyFloat64Writer out(path,shape);
    for(std::size_t offset=0; offset<buffer.element_count();)
    {
        const auto n=(std::min)(staging.count,buffer.element_count()-offset);
        // Synchronous copy: the reusable pinned buffer never outlives a transfer,
        // including on a write/close exception. No second whole-table allocation.
        c.make_current();
        RAINBOW_CUDA_CHECK(cuMemcpyDtoH(staging.data,buffer.address()+offset*sizeof(double),n*sizeof(double)));
        out.append(std::span<const double>{staging.data,n});
        offset+=n;
    }
    out.finish();
}
}
PhaseCdf::PhaseCdf(const CudaContext& c)
    : context_(c),module_(c),edges_(c),phi_(c),theta_(c),column_sums_(c),total_(c),
      partials_(c),input_summary_(c),audit_summary_(c),
      fine_edges_(c),fine_values_(c),filtered_values_(c),coarse_values_(c),
      moment_partials_(c),moment_result_(c) {}
PhaseCdf::~PhaseCdf() noexcept
{
    if(pending_)
    {
        static_cast<void>(detail::report_cuda_cleanup_result(cuCtxSetCurrent(context_.handle()),"CDF cleanup context"));
        static_cast<void>(detail::report_cuda_cleanup_result(cuStreamSynchronize(context_.stream()),"CDF cleanup sync"));
    }
}
void PhaseCdf::load_module(const std::filesystem::path& path)
{
    module_.load_fatbin(path);
    static_cast<void>(module_.find_function("rainbow_phase_cdf_abi_v1"));
    static_cast<void>(module_.find_function("rainbow_phase_storage_abi_v2"));
    // Diagnose a stale CDF module at load time, not after an expensive trace.
    static_cast<void>(module_.find_function("phase_check_density_failures_only"));
}
void PhaseCdf::synchronize()
{
    context_.make_current();
    RAINBOW_CUDA_CHECK(cuStreamSynchronize(context_.stream()));
    pending_=false;
}
void PhaseCdf::launch(const char* name, void** args, unsigned gx,unsigned gy,unsigned tx,unsigned ty)
{
    context_.make_current();
    const auto f=module_.find_function(name);
    pending_=true;
    RAINBOW_CUDA_CHECK(cuLaunchKernel(f,gx,gy,1,tx,ty,1,0,context_.stream(),args,nullptr));
}
void PhaseCdf::allocate_partials(std::uint64_t count)
{
    partial_count_=static_cast<std::uint32_t>((std::min)(std::uint64_t{4096},(count+255)/256));
    partials_.close(); partials_.allocate(partial_count_);
    if(input_summary_.is_empty()) input_summary_.allocate(1);
    if(audit_summary_.is_empty()) audit_summary_.allocate(1);
}
void PhaseCdf::reduce_to(DeviceBuffer<PhaseCdfReduction>& output)
{
    const auto* in=reinterpret_cast<const PhaseCdfReduction*>(partials_.address());
    auto* out=output.data();
    void* args[]={&in,&partial_count_,&out};
    launch("phase_reduce_summary",args,1);
}
void PhaseCdf::validate_trace(const RaindropTracer& trace)
{
    if(!module_.is_loaded() || !trace.has_result() || !trace.is_unpolarized() ||
       trace.vertices().element_count()!=trace.second_input_fields().element_count() ||
       trace.vertices().is_empty())
        throw std::invalid_argument("Trace validation requires a completed two-input trace and CDF module.");
    synchronize(); allocate_partials(trace.vertices().element_count());
    TraceValidationParams p{};
    p.vertices=reinterpret_cast<const OutgoingVertex*>(trace.vertices().address());
    p.second=reinterpret_cast<const Field32*>(trace.second_input_fields().address());
    p.count=trace.vertices().element_count();p.partials=partials_.data();
    void* args[]={&p};launch("phase_validate_trace",args,partial_count_);
    reduce_to(input_summary_);synchronize();
    PhaseCdfReduction status{};input_summary_.download({&status,1});
    if(status.incomplete) throw std::runtime_error("Trace contains failed vertices: "+describe(status));
}
void PhaseCdf::build(const PatchOptics& optics,std::uint32_t nt,std::uint32_t np,
                     PhaseDensityStage stage,const PhaseCdfPolicy& policy,const PhaseStorageSettings& storage)
{
    valid_=false;
    if(!optics.has_result() || !optics.is_unpolarized() || stage==PhaseDensityStage::Scalar)
        throw std::invalid_argument("Require an unpolarized optical result and optical stage.");
    const auto* grid=optics.direction_grid();
    if(!grid || grid->theta_count!=nt || grid->phi_count!=np)
        throw std::invalid_argument("Optical angular-grid metadata does not match the CDF dimensions.");
    if((stage==PhaseDensityStage::Focal || stage==PhaseDensityStage::Diffraction) && !optics.has_wave_result())
        throw std::invalid_argument("Wave result is absent.");
    PhaseDensityView v{};v.stage=stage;v.count=optics.results().element_count();
    v.optical=reinterpret_cast<const PatchOpticalResult*>(optics.results().address());
    v.focal=reinterpret_cast<const FocalOpticalResult*>(optics.focal_results().address());
    v.diffraction=reinterpret_cast<const DiffractionResult*>(optics.diffraction_results().address());
    if(v.focal && optics.focal_results().element_count()!=v.count)
        throw std::invalid_argument("Focal array length mismatch.");
    if(v.diffraction && optics.diffraction_results().element_count()!=v.count)
        throw std::invalid_argument("Diffraction array length mismatch.");
    build(v,nt,np,policy,storage);
}
void PhaseCdf::build(PhaseDensityView input,std::uint32_t nt,std::uint32_t np,const PhaseCdfPolicy& policy,const PhaseStorageSettings& storage)
{
    valid_=false;validate_policy(policy);
    const auto out_nt=storage.theta_count?storage.theta_count:nt;
    const auto out_np=storage.phi_count?storage.phi_count:np;
    if(!out_nt || !out_np || !nt || !np || out_nt>nt || out_np>np || nt%out_nt || np%out_np)
        throw std::invalid_argument("CDF storage grid must be an integer coarsening of the query grid.");
    const double sigma=storage.gaussian_sigma_degrees*phase_cdf_pi/180.0;
    if(!(std::isfinite(sigma) && sigma>=0 && std::isfinite(storage.gaussian_support_sigma)
        && storage.gaussian_support_sigma>=2 && storage.gaussian_support_sigma<=8
        && sigma*storage.gaussian_support_sigma<phase_cdf_pi
        && std::isfinite(storage.maximum_coarsening_tv) && storage.maximum_coarsening_tv>=0
        && storage.maximum_coarsening_tv<=1))
        throw std::invalid_argument("Invalid storage Gaussian / coarsening error policy.");
    if(sigma>0 && (nt>32767 || np>32767))
        throw std::invalid_argument("Storage Gaussian uses the solver's <=32767 angular-index range.");
    storage_statistics_={};
    storage_statistics_.gaussian_underresolved=sigma>0 &&
        ((phase_cdf_pi/double(nt)>0.5*sigma)||(2*phase_cdf_pi/double(np)>0.5*sigma));
    if(storage_statistics_.gaussian_underresolved && !policy.allow_underresolved)
        throw std::invalid_argument("Storage Gaussian is underresolved on the QUERY grid; use a finer query grid or explicit --allow-underresolved.");
    const auto count=std::uint64_t(nt)*np;
    // Same angular-index range as the existing solver. Tiles must fit gridDim.y.
    if(!module_.is_loaded() || nt==0 || np==0 || count>0xffffffffull ||
       nt>32u*65535u || np>0x7fffffffu || input.count!=count ||
       static_cast<unsigned>(input.stage)>static_cast<unsigned>(PhaseDensityStage::Diffraction))
        throw std::invalid_argument("Invalid CDF dimensions/module/stage.");
    if((input.stage==PhaseDensityStage::Scalar && !input.scalar) ||
       (input.stage!=PhaseDensityStage::Scalar && !input.optical) ||
       ((input.stage==PhaseDensityStage::Focal || input.stage==PhaseDensityStage::Diffraction) && !input.focal) ||
       (input.stage==PhaseDensityStage::Diffraction && !input.diffraction))
        throw std::invalid_argument("Missing source buffer.");
    synchronize();nt_=nt;np_=np;
    edges_.close();phi_.close();theta_.close();column_sums_.close();total_.close();
    fine_edges_.close();fine_values_.close();filtered_values_.close();coarse_values_.close();
    moment_partials_.close();moment_result_.close();
    allocate_partials(count);
    moment_partials_.allocate(partial_count_);moment_result_.allocate(1);
    host_edges_=make_phase_u_edges(nt);
    edges_.allocate(host_edges_.size());
    // This O(Ntheta) coordinate metadata is not an optical readback.
    upload_host_to_device_sync(context_, edges_.address(), host_edges_.data(), edges_.byte_size());

    PhaseCdfBuildParams p{};
    p.input=input;p.u_edges=edges_.data();p.input_summary=input_summary_.data();
    p.partials=partials_.data();p.output_summary=audit_summary_.data();
    p.theta_count=nt;p.phi_count=np;
    void* args[]={&p};
    launch("phase_check_density",args,partial_count_);reduce_to(input_summary_);synchronize();
    input_summary_.download({&input_statistics_,1});
    if(input_statistics_.invalid_values || input_statistics_.incomplete ||
       !(input_statistics_.maximum>0) ||
       (!policy.allow_underresolved && input_statistics_.underresolved))
    {
        std::ostringstream failure;
        failure.imbue(std::locale::classic());
        failure << std::setprecision(17)
            << "CDF input rejected: " << describe(input_statistics_) << '\n'
            << "[cdf-input] stage=" << stage_name(input.stage)
            << " directions=" << count << " allow_underresolved=" << policy.allow_underresolved
            << " CDF_accumulation_started=0\n"
            << "[cdf-input] original first_problem includes underresolution warnings; "
               "first_hard_problem below does not.\n";
        // Extra scans run ONLY when the original gate rejects the input.
        // Reuse the already allocated reduction scratch. No ray trace, patch
        // query, optical evaluation, diffraction or bulk readback is repeated.
        try
        {
            const unsigned first=input.stage==PhaseDensityStage::Scalar ? 0u : 1u;
            const unsigned last=static_cast<unsigned>(input.stage);
            for(unsigned s=first;s<=last;++s)
            {
                auto inspect=p;
                inspect.input.stage=static_cast<PhaseDensityStage>(s);
                void* inspect_args[]={&inspect};
                launch("phase_check_density_failures_only",inspect_args,partial_count_);
                reduce_to(audit_summary_);
                synchronize();
                PhaseCdfReduction status{};
                audit_summary_.download({&status,1});
                failure << "[cdf-input] stage=" << stage_name(inspect.input.stage)
                    << " invalid_values=" << status.invalid_values
                    << " incomplete=" << status.incomplete
                    << " underresolved=" << status.underresolved
                    << " first_hard_problem=";
                if(status.first_problem==~std::uint64_t{0}) failure << "none";
                else failure << status.first_problem;
                failure << '\n';
                if(status.first_problem!=~std::uint64_t{0})
                    append_failure_record(failure,context_,inspect.input,
                                          status.first_problem,nt,np);
            }
        }
        catch(const std::exception& e)
        {
            // The original input rejection is still fatal. Do not replace it
            // by a successful return if the optional inspection itself fails.
            failure << "[cdf-input] additional_inspection_failed=" << e.what() << '\n';
        }
        throw std::runtime_error(failure.str());
    }
    // The original optical gate ABOVE is mandatory: smoothing must never hide
    // an unavailable contribution. The source statistics keep their old meaning.
    PhaseStorageParams sp{};sp.input=input;sp.nt=nt;sp.np=np;
    sp.fine_edges=edges_.data();sp.divisor=input_statistics_.maximum;
    sp.out_nt=out_nt;sp.out_np=out_np;
    const auto source=measure(sp,"phase_storage_measure");
    if(source.weight_underflow || !(std::isfinite(source.mass)&&source.mass>0&&std::isfinite(source.axial)))
        throw std::runtime_error("Invalid source angular moment." + describe_moment(source,sp)
                + inspect_moment_edges(context_,sp,host_edges_));
    storage_statistics_.g_source=source.axial/source.mass;
    storage_statistics_.g_filtered=storage_statistics_.g_source;
    const bool transformed=sigma>0 || out_nt!=nt || out_np!=np;
    if(transformed)
    {
        fine_edges_.allocate(host_edges_.size());
        upload_host_to_device_sync(context_, fine_edges_.address(), host_edges_.data(), fine_edges_.byte_size());
        sp.fine_edges=fine_edges_.data();
        fine_values_.allocate(count);sp.output=fine_values_.data();
        void* storage_args[]={&sp};
        launch("phase_storage_extract",storage_args,partial_count_);
        sp.fine_values=fine_values_.data();
        if(sigma>0)
        {
            filtered_values_.allocate(count);sp.output=filtered_values_.data();
            sp.sigma_rad=sigma;sp.support_sigma=storage.gaussian_support_sigma;
            launch("phase_storage_gaussian",storage_args,partial_count_);
            sp.fine_values=filtered_values_.data();
        }
        sp.input={};sp.input.scalar=sp.fine_values;sp.input.count=count;sp.divisor=1;
        const auto filtered=measure(sp,"phase_storage_measure");
        if(filtered.weight_underflow || !(std::isfinite(filtered.mass)&&filtered.mass>0&&std::isfinite(filtered.axial)))
            throw std::runtime_error("Storage filtering produced an invalid angular moment." + describe_moment(filtered,sp)
                + inspect_moment_edges(context_,sp,host_edges_));
        storage_statistics_.g_filtered=filtered.axial/filtered.mass;
        storage_statistics_.gaussian_integral_relative_change=filtered.mass/source.mass-1.0;
        if(out_nt!=nt || out_np!=np)
        {
            coarse_values_.allocate(std::size_t(out_nt)*out_np);sp.output=coarse_values_.data();
            launch("phase_storage_coarsen",storage_args,partial_count_);
            std::vector<double> coarse_edges(std::size_t(out_nt)+1);
            for(std::uint32_t i=0;i<=out_nt;++i)coarse_edges[i]=host_edges_[i*(nt/out_nt)];
            host_edges_=std::move(coarse_edges);
            // Previous kernels only read fine_edges_, not edges_.
            synchronize();edges_.close();edges_.allocate(host_edges_.size());
            upload_host_to_device_sync(context_, edges_.address(), host_edges_.data(), edges_.byte_size());
            auto cp=sp;cp.input={};cp.input.scalar=coarse_values_.data();cp.input.count=std::uint64_t(out_nt)*out_np;
            cp.nt=out_nt;cp.np=out_np;cp.fine_edges=edges_.data();
            const auto coarse=measure(cp,"phase_storage_measure");
            if(coarse.weight_underflow || !(std::isfinite(coarse.mass)&&coarse.mass>0&&std::isfinite(coarse.axial)))
                throw std::runtime_error("Storage aggregation produced an invalid angular moment." + describe_moment(coarse,cp)
                + inspect_moment_edges(context_,cp,host_edges_));
            storage_statistics_.aggregation_integral_relative_change=coarse.mass/filtered.mass-1.0;
            if(::fabs(storage_statistics_.aggregation_integral_relative_change)>1e-10)
                throw std::runtime_error("Storage aggregation failed mass-conservation check.");
            sp.coarse_values=coarse_values_.data();sp.fine_mass=filtered.mass;sp.coarse_mass=coarse.mass;
            storage_statistics_.coarsening_tv=measure(sp,"phase_storage_measure_tv").tv;
            if(!std::isfinite(storage_statistics_.coarsening_tv) ||
                storage_statistics_.coarsening_tv>storage.maximum_coarsening_tv)
                throw std::runtime_error("Storage coarsening exceeds configured TV limit: "+std::to_string(storage_statistics_.coarsening_tv));
            input=cp.input;
        }
        else input=sp.input;
        nt=out_nt;np=out_np;nt_=nt;np_=np;
        p.input=input;p.u_edges=edges_.data();p.theta_count=nt;p.phi_count=np;
        launch("phase_check_density",args,partial_count_);reduce_to(input_summary_);synchronize();
        PhaseCdfReduction processed{};input_summary_.download({&processed,1});
        if(processed.invalid_values || processed.incomplete || !(processed.maximum>0))
            throw std::runtime_error("Invalid processed density: "+describe(processed));
    }
    phi_.allocate(std::size_t(np)+1);theta_.allocate(std::size_t(np)*(std::size_t(nt)+1));
    column_sums_.allocate(np);total_.allocate(1);
    p.phi_cdf=phi_.data();p.theta_cdf=theta_.data();p.column_sums=column_sums_.data();p.total=total_.data();
    launch("phase_make_weights",args,(np+31u)/32u,(nt+31u)/32u,32,8);
    launch("phase_scan_conditionals",args,(np+255u)/256u);
    launch("phase_scan_marginal",args,1,1,1,1);
    launch("phase_normalize",args,partial_count_);
    launch("phase_audit_cdf",args,partial_count_);reduce_to(audit_summary_);synchronize();
    audit_summary_.download({&audit_statistics_,1});
    if(audit_statistics_.malformed || audit_statistics_.weight_underflow ||
       !std::isfinite(audit_statistics_.l1_error) ||
       !std::isfinite(audit_statistics_.lost_mass) ||
       audit_statistics_.l1_error>policy.maximum_l1_error ||
       audit_statistics_.lost_mass>policy.maximum_lost_mass)
        throw std::runtime_error("CDF rounding/structure audit failed: "+describe(audit_statistics_));
    sp.phi_cdf=phi_.data();sp.theta_cdf=theta_.data();sp.coarse_edges=edges_.data();
    sp.out_nt=nt;sp.out_np=np;
    const auto saved=measure(sp,"phase_storage_measure_cdf");
    if(!(std::isfinite(saved.mass)&&saved.mass>0&&std::isfinite(saved.axial)) || ::fabs(saved.mass-1)>1e-10)
        throw std::runtime_error("Invalid probability mass / moment in stored CDF.");
    storage_statistics_.cdf_mass=saved.mass;
    storage_statistics_.g_stored=saved.axial/saved.mass;
    if(!(storage_statistics_.g_stored>-1 && storage_statistics_.g_stored<1))
        throw std::runtime_error("Stored first moment is outside the nondegenerate HG range.");
    // All launches are synchronized by measure(). Only CDFs survive the build.
    fine_values_.close();filtered_values_.close();coarse_values_.close();fine_edges_.close();
    valid_=true;
}
PhaseMomentSum PhaseCdf::measure(PhaseStorageParams p,const char* kernel)
{
    p.moment_partials=moment_partials_.data();
    void* args[]={&p};launch(kernel,args,partial_count_);
    auto* in=moment_partials_.data();auto* out=moment_result_.data();
    void* reduction[]={&in,&partial_count_,&out};
    launch("phase_storage_reduce_moments",reduction,1);synchronize();
    PhaseMomentSum result{};moment_result_.download({&result,1});return result;
}
void PhaseCdf::write_arrays(const std::filesystem::path& directory) const
{
    if(!valid_) throw std::logic_error("No validated CDF to save.");
    context_.make_current();
    PinnedDoubles staging(1u<<20); // 8 MiB, reused for every output array.
    const std::array<std::uint64_t,1> a{std::uint64_t(np_)+1};
    const std::array<std::uint64_t,2> b{np_,std::uint64_t(nt_)+1};
    const std::array<std::uint64_t,1> e{std::uint64_t(nt_)+1};
    write_device_npy(context_,phi_,directory/"phi_cdf.npy",a,staging);
    write_device_npy(context_,theta_,directory/"theta_given_phi_cdf.npy",b,staging);
    NpyFloat64Writer out(directory/"u_edges.npy",e);out.append(host_edges_);out.finish();
}
} // namespace rainbow
