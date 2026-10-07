#include <rainbow/phase_cdf.hpp>
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
      partials_(c),input_summary_(c),audit_summary_(c) {}
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
                     PhaseDensityStage stage,const PhaseCdfPolicy& policy)
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
    build(v,nt,np,policy);
}
void PhaseCdf::build(PhaseDensityView input,std::uint32_t nt,std::uint32_t np,const PhaseCdfPolicy& policy)
{
    valid_=false;validate_policy(policy);
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
    allocate_partials(count);
    host_edges_=make_phase_u_edges(nt);
    edges_.allocate(host_edges_.size());
    // This O(Ntheta) coordinate metadata is not an optical readback.
    RAINBOW_CUDA_CHECK(cuMemcpyHtoD(edges_.address(),host_edges_.data(),edges_.byte_size()));
    phi_.allocate(std::size_t(np)+1);theta_.allocate(std::size_t(np)*(std::size_t(nt)+1));
    column_sums_.allocate(np);total_.allocate(1);
    PhaseCdfBuildParams p{};
    p.input=input;p.u_edges=edges_.data();p.input_summary=input_summary_.data();
    p.phi_cdf=phi_.data();p.theta_cdf=theta_.data();p.column_sums=column_sums_.data();
    p.total=total_.data();p.partials=partials_.data();p.output_summary=audit_summary_.data();
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
    valid_=true;
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
