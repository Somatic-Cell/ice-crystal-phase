#include <ice_crystal/phase_generate.hpp>
#include <rainbow/npy_writer.hpp>
#if defined(ICE_PHASE_HAS_CUDA)
#include <ice_crystal/phase_cuda.hpp>
#include <rainbow/phase_cdf.hpp>
#endif
#include <algorithm>
#include <charconv>
#include <iostream>
#include <memory>
#include <set>
#include <string_view>
#include <limits>
#include <stdexcept>

namespace
{
using namespace iceCrystal;
std::filesystem::path path_utf8(std::string_view s)
{std::u8string v;for(unsigned char c:s)v+=static_cast<char8_t>(c);return std::filesystem::path(v);}
template<class T>T number(std::string_view s)
{
    T v{};const auto r=std::from_chars(s.data(),s.data()+s.size(),v);
    if(r.ec!=std::errc{}||r.ptr!=s.data()+s.size())throw std::invalid_argument("Malformed numeric option.");
    return v;
}
struct Options
{
    PhaseGenerationSettings settings;
    std::filesystem::path out,orientations,modules;
    int device=0;
};
Options parse(int argc,char** argv)
{
    Options o;auto& s=o.settings;Vec3 ki{};std::uint32_t nt=0,np=0;
    o.modules=std::filesystem::absolute(path_utf8(argv[0])).parent_path()/"modules";
    std::set<std::string> seen;
    for(int i=1;i<argc;++i)
    {
        const std::string key=argv[i];if(!seen.insert(key).second)throw std::invalid_argument("Duplicate option: "+key);
        auto next=[&]()->std::string_view{if(++i>=argc)throw std::invalid_argument("Missing value: "+key);return argv[i];};
        if(key=="--out")o.out=path_utf8(next());
        else if(key=="--orientations")o.orientations=path_utf8(next());
        else if(key=="--modules")o.modules=path_utf8(next());
        else if(key=="--radius-mm")s.prism.circumradius_mm=number<double>(next());
        else if(key=="--length-mm")s.prism.length_mm=number<double>(next());
        else if(key=="--ki"){ki.x=number<double>(next());ki.y=number<double>(next());ki.z=number<double>(next());}
        else if(key=="--ior")s.interior_index=number<double>(next());
        else if(key=="--external-ior")s.exterior_index=number<double>(next());
        else if(key=="--wavelength-nm")s.wavelength_nm=number<double>(next());
        else if(key=="--theta")nt=number<std::uint32_t>(next());
        else if(key=="--phi")np=number<std::uint32_t>(next());
        else if(key=="--cdf-theta")s.out_nt=number<std::uint32_t>(next());
        else if(key=="--cdf-phi")s.out_np=number<std::uint32_t>(next());
        else if(key=="--batch-size")s.batch_size=number<std::size_t>(next());
        else if(key=="--max-output-records")s.max_output_records=number<std::size_t>(next());
        else if(key=="--max-internal-hits")s.max_internal_hits=number<std::uint32_t>(next());
        else if(key=="--tail-tolerance")s.tail_tolerance=number<double>(next());
        else if(key=="--max-unresolved-fraction")s.max_unresolved_fraction=number<double>(next());
        else if(key=="--max-balance-error")s.max_balance_error=number<double>(next());
        else if(key=="--max-histogram-relative-error")s.max_histogram_relative_error=number<double>(next());
        else if(key=="--cdf-max-l1")s.cdf_policy.max_l1_error=number<double>(next());
        else if(key=="--cdf-max-lost-mass")s.cdf_policy.max_lost_mass=number<double>(next());
        else if(key=="--max-coarsening-tv")s.cdf_policy.max_coarsening_tv=number<double>(next());
        else if(key=="--device")o.device=number<int>(next());
        else throw std::invalid_argument("Unknown option: "+key);
    }
    for(const char* k:{"--out","--orientations","--radius-mm","--length-mm","--ki","--ior","--wavelength-nm","--theta","--phi"})
        if(!seen.contains(k))throw std::invalid_argument(std::string("Required option: ")+k);
#if !defined(ICE_PHASE_HAS_CUDA)
    if(seen.contains("--device")||seen.contains("--modules"))throw std::invalid_argument("CPU reference binary does not accept CUDA settings.");
#endif
    if(o.device<0)throw std::invalid_argument("Negative CUDA device ordinal.");
    s.grid=make_phase_grid(nt,np,ki);
    if(!seen.contains("--cdf-theta"))s.out_nt=nt;
    if(!seen.contains("--cdf-phi"))s.out_np=np;
    validate_generation_settings(s);return o;
}
int run(const Options& o)
{
    const auto& s=o.settings;const auto plan=read_orientation_plan(o.orientations);CompositionAudit audit;
    // Preflight all nodes before producing any output.
    for(const auto& node:plan.nodes)
    {
        const double area=node_area(s,node),weight=node.weight*area;
        if(!(weight>0)||!finite_value(weight)||!(weight/double(node.samples)>0))throw std::runtime_error("Orientation area weight over/underflow.");
        audit.expected_area.add(weight);
    }
    rainbow::DatasetDirectory directory(o.out);
    try
    {
#if defined(ICE_PHASE_HAS_CUDA)
        rainbow::CudaContext context(o.device);HexPrismTracer tracer(context);
        tracer.load_module(o.modules/"hex_trace.fatbin");
        PhaseMassAccumulatorCuda histogram(context,s.grid);histogram.load_module(o.modules/"phase_accumulate.fatbin");
#else
        PhaseMassAccumulator histogram(s.grid);
#endif
        for(std::size_t pose=0;pose<plan.nodes.size();++pose)
        {
            const auto& node=plan.nodes[pose];const auto c=node_trace_settings(s,node);
            const double scale=node.weight*node_area(s,node)/double(node.samples);
            for(std::uint64_t first=0;first<node.samples;)
            {
                const auto count=static_cast<std::size_t>(std::min<std::uint64_t>(s.batch_size,node.samples-first));
#if defined(ICE_PHASE_HAS_CUDA)
                tracer.trace(c,first,count,s.max_output_records);
                const auto rays=tracer.download_audits();
                record_batch_audit(audit,summarize_phase_rays(rays,s.max_balance_error),scale);
                histogram.add(tracer.outgoing(),node.rotation,scale);
#else
                const auto batch=trace_batch_cpu(c,first,count,s.max_output_records);
                record_batch_audit(audit,summarize_phase_rays(batch.audits,s.max_balance_error),scale);
                histogram.add(batch.outgoing,node.rotation,scale);
#endif
                first+=count;
            }
            ++audit.orientations_done;
            std::cout<<"orientation="<<pose+1<<'/'<<plan.nodes.size()<<" rays="<<audit.trace.rays<<'\n';
        }
        CdfReport report;
#if defined(ICE_PHASE_HAS_CUDA)
        histogram.finish();audit.histogram=histogram.statistics();audit.grid_moments=histogram.moments();
        validate_composition(audit,s,plan);
        rainbow::PhaseCdf cdf(context);cdf.load_module(o.modules/"phase_cdf.fatbin");
        rainbow::PhaseDensityView input{};input.stage=rainbow::PhaseDensityStage::Scalar;
        input.count=histogram.density().element_count();input.scalar=reinterpret_cast<const double*>(histogram.density().address());
        rainbow::PhaseCdfPolicy policy{};policy.maximum_l1_error=s.cdf_policy.max_l1_error;policy.maximum_lost_mass=s.cdf_policy.max_lost_mass;
        rainbow::PhaseStorageSettings storage{};storage.theta_count=s.out_nt;storage.phi_count=s.out_np;
        storage.gaussian_sigma_degrees=0;storage.maximum_coarsening_tv=s.cdf_policy.max_coarsening_tv;
        cdf.build(input,s.grid.nt,s.grid.np,policy,storage);
        const auto& p=cdf.storage_statistics();const auto& q=cdf.audit_statistics();
        report={p.g_source,p.g_stored,p.coarsening_tv,p.cdf_mass,q.l1_error,q.lost_mass,q.maximum_cell_error,p.aggregation_integral_relative_change,q.lost_cells};
        cdf.write_arrays(directory.staging());const char* backend="cuda_trace_atomic_histogram_existing_rainbow_PhaseCdf";
#else
        const auto mass=histogram.masses();audit.histogram=histogram.statistics();audit.grid_moments=measure_mass(s.grid,mass);
        validate_composition(audit,s,plan);
        const auto cdf=build_host_cdf(s.grid,mass,s.out_nt,s.out_np,s.cdf_policy);report=cdf.report;
        cdf.write_arrays(directory.staging());const char* backend="cpu_reference_compensated_histogram_host_mass_CDF";
#endif
        write_phase_metadata(directory.staging()/"metadata.json",s,plan,audit,report,backend);
        directory.commit();
        std::cout.precision(17);
        std::cout<<"saved="<<o.out<<" g="<<report.g_stored<<" Csca_mm2="<<audit.expected_area.value
                 <<" unresolved_fraction="<<audit.unresolved_area.value/audit.expected_area.value
                 <<" histogram_relative_error="<<audit.histogram_relative_error<<" cdf_l1="<<report.l1_error<<'\n';
        return 0;
    }
    catch(const std::exception& e)
    {
        write_failure(directory.staging()/"failure.json",e.what(),audit);
        std::cerr<<"Generation failed; retained "<<directory.staging()<<": "<<e.what()<<'\n';return 2;
    }
}
int run_main(int argc,char** argv)
{
    try
    {
        if(argc==2&&std::string_view(argv[1])=="--help")
        {
            std::cout<<"ice_generate_phase[_cpu] --out DIR --orientations CSV\n"
                " --radius-mm R --length-mm L --ki X Y Z --ior N --wavelength-nm NM\n"
                " --theta NT --phi NP [--cdf-theta CT --cdf-phi CP]\n"
                " [--batch-size N --max-output-records N --max-internal-hits N --tail-tolerance E]\n"
                " [--max-unresolved-fraction E --max-balance-error E --max-histogram-relative-error E]\n"
                " [--cdf-max-l1 E --cdf-max-lost-mass E --max-coarsening-tv E]\n"
                " CUDA binary only: [--modules DIR --device N]\n"
                "No implicit material dispersion, smoothing, symmetry folding, or CPU fallback.\n";return 0;
        }
        return run(parse(argc,argv));
    }
    catch(const std::exception& e){std::cerr<<"ice_generate_phase: "<<e.what()<<'\n';return 1;}
}
}
#if defined(_WIN32)
int wmain(int argc,wchar_t** argv)
{
    try
    {
        std::vector<std::string> strings;strings.reserve(static_cast<std::size_t>(argc));
        for(int i=0;i<argc;++i){const auto u=std::filesystem::path(argv[i]).u8string();strings.emplace_back(reinterpret_cast<const char*>(u.data()),u.size());}
        std::vector<char*> args;for(auto& s:strings)args.push_back(s.data());return run_main(argc,args.data());
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
#else
int main(int argc,char** argv){return run_main(argc,argv);}
#endif
