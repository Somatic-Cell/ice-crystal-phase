#include <ice_crystal/atlas.hpp>
#include <ice_crystal/ice_optics.hpp>
#include <rainbow/npy_writer.hpp>
#if defined(ICE_ATLAS_HAS_CUDA)
#include <ice_crystal/atlas_cuda.hpp>
#endif
#include <array>
#include <charconv>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <set>
#include <stdexcept>
#include <string>
#ifndef ICE_ATLAS_SOURCE_FINGERPRINT
#define ICE_ATLAS_SOURCE_FINGERPRINT "unknown"
#endif
namespace
{
using namespace iceCrystal;
struct Options
{
    atlas::Settings s;
    std::filesystem::path out,modules;
    int device=0;
    bool ice_material=true;
    optics::IceState optical_state;
};
double real(std::string_view v)
{
    double x=0;const auto r=std::from_chars(v.data(),v.data()+v.size(),x);
    if(r.ec!=std::errc{}||r.ptr!=v.data()+v.size()||!finite_value(x))throw std::invalid_argument("Invalid finite number.");
    return x;
}
std::uint32_t integer(std::string_view v)
{
    std::uint32_t x=0;const auto r=std::from_chars(v.data(),v.data()+v.size(),x);
    if(r.ec!=std::errc{}||r.ptr!=v.data()+v.size())throw std::invalid_argument("Invalid integer.");
    return x;
}
std::filesystem::path utf8_path(std::string_view s)
{
    return std::filesystem::path(std::u8string_view(reinterpret_cast<const char8_t*>(s.data()),s.size()));
}
Options parse(int argc,char** argv)
{
    Options o;o.modules=std::filesystem::absolute(utf8_path(argv[0])).parent_path()/"modules";
    std::set<std::string> seen;
    for(int i=1;i<argc;++i)
    {
        const std::string key=argv[i];
        if(!seen.insert(key).second||i+1>=argc)throw std::invalid_argument("Duplicate option or missing argument.");
        const std::string value=argv[++i];auto& s=o.s;
        if(key=="--out")o.out=utf8_path(value);
        else if(key=="--modules")o.modules=utf8_path(value);
        else if(key=="--device") {const auto v=integer(value);if(v>2147483647u)throw std::invalid_argument("Device range.");o.device=int(v);}
        else if(key=="--radius-mm")s.prism.circumradius_mm=real(value);
        else if(key=="--length-mm")s.prism.length_mm=real(value);
        else if(key=="--material")
        {
            if(value=="ice")o.ice_material=true;
            else if(value=="constant")o.ice_material=false;
            else throw std::invalid_argument("Material must be ice or constant.");
        }
        else if(key=="--ior")s.interior_index=real(value);
        else if(key=="--exterior-ior")s.exterior_index=real(value);
        else if(key=="--wavelength-nm")s.wavelength_nm=real(value);
        else if(key=="--incident-theta-intervals")s.incident_theta_intervals=integer(value);
        else if(key=="--incident-phi-intervals")s.incident_phi_intervals=integer(value);
        else if(key=="--theta")s.theta_count=integer(value);
        else if(key=="--phi")s.phi_count=integer(value);
        else if(key=="--area-subdivisions")s.area_subdivisions=integer(value);
        else if(key=="--batch-rays")s.batch_rays=integer(value);
        else if(key=="--max-work-mib")s.max_work_bytes=std::size_t(integer(value))*1024*1024;
        else if(key=="--max-internal-hits")s.max_internal_hits=integer(value);
        else if(key=="--tail-tolerance")s.tail_tolerance=real(value);
        else if(key=="--max-unresolved-fraction")s.max_unresolved_fraction=real(value);
        else if(key=="--max-balance-error")s.max_balance_error=real(value);
        else throw std::invalid_argument("Unknown option: "+key);
    }
#if !defined(ICE_ATLAS_HAS_CUDA)
    if(seen.count("--device")||seen.count("--modules"))throw std::invalid_argument("CPU binary does not accept CUDA settings.");
#endif
    if(o.out.empty())throw std::invalid_argument("--out required.");
    if(!seen.count("--wavelength-nm"))throw std::invalid_argument("--wavelength-nm is required (vacuum nm).");
    if(o.ice_material)
    {
        if(seen.count("--ior"))throw std::invalid_argument("Ice index is evaluated from wavelength; --ior requires explicit --material constant.");
        o.optical_state=optics::evaluate_ice_wb2008(o.s.wavelength_nm);
        o.s.interior_index=o.optical_state.n;
    }
    else
    {
        if(!seen.count("--ior"))throw std::invalid_argument("Constant material requires --ior.");
        o.optical_state={o.s.wavelength_nm,o.s.interior_index,0,0,0,0};
    }
    atlas::validate(o.s);return o;
}
void vec(std::ostream& o,Vec3 v){o<<'['<<v.x<<','<<v.y<<','<<v.z<<']';}
int run(int argc,char** argv)
{
    if(argc==2&&std::string_view(argv[1])=="--help")
    {
        std::cout<<"ice_precompute_atlas[_cpu] --out DIR --radius-mm R --length-mm L --wavelength-nm NM\n"
          " --material ice (default, WB2008 n(lambda), 200..1100nm) | --material constant --ior N\n"
          " --incident-theta-intervals T --incident-phi-intervals P --theta NT --phi NP --area-subdivisions S\n"
          " [--max-internal-hits 2048 --tail-tolerance 1e-12 --batch-rays 64 --max-work-mib 512]\n"
          " [--max-unresolved-fraction 1e-10 --max-balance-error 1e-10 --exterior-ior 1]\n"
          "CUDA only: [--modules DIR --device 0]. No random seeds, smoothing, or CPU fallback.\n";return 0;
    }
    const auto o=parse(argc,argv);const auto& s=o.s;
    std::cout<<std::setprecision(17)<<"material="<<(o.ice_material?"ice_wb2008_real_only":"constant")
             <<" wavelength_nm="<<s.wavelength_nm<<" n="<<s.interior_index
             <<" kappa_evaluated_not_used="<<o.optical_state.kappa<<'\n';
    rainbow::DatasetDirectory directory(o.out);
    try
    {
#if defined(ICE_ATLAS_HAS_CUDA)
        rainbow::CudaContext context(o.device);atlas::CudaTracer tracer(context,o.modules/"atlas_trace.fatbin");
        constexpr const char* backend="cuda_ordered_per_ray_histogram";
#else
        constexpr const char* backend="cpu_ordered_per_ray_histogram";
#endif
        const std::uint64_t ni=std::uint64_t(s.incident_theta_intervals)+1,nj=std::uint64_t(s.incident_phi_intervals)+1;
        const std::array<std::uint64_t,4> shape{ni,nj,s.theta_count,s.phi_count};
        rainbow::NpyFloat64Writer mass(directory.staging()/"mass_mm2.npy",shape);
        std::ofstream meta(directory.staging()/"metadata.json",std::ios::binary|std::ios::trunc);
        meta.exceptions(std::ios::badbit|std::ios::failbit);meta.imbue(std::locale::classic());meta<<std::setprecision(17);
        meta<<"{\n\"schema\":\"ice.single_particle_atlas.v1\",\"complete\":true,\"backend\":\""<<backend<<"\",\n"
          <<"\"source_fingerprint_at_configure\":\""<<ICE_ATLAS_SOURCE_FINGERPRINT<<"\",\n"
          <<"\"dtype\":\"<f8\",\"order\":\"C\",\"mass_axes\":[\"incident_theta\",\"incident_phi\",\"theta\",\"phi\"],\n"
          <<"\"incident_theta_intervals\":"<<s.incident_theta_intervals<<",\"incident_phi_intervals\":"<<s.incident_phi_intervals<<",\n"
          <<"\"theta_count\":"<<s.theta_count<<",\"phi_count\":"<<s.phi_count<<",\n"
          <<"\"geometry\":{\"shape\":\"regular_hexagonal_prism\",\"circumradius_mm\":"<<s.prism.circumradius_mm
          <<",\"full_length_mm\":"<<s.prism.length_mm<<",\"body_long_axis\":[0,1,0]},\n"
          <<"\"material\":";
        optics::write_material_json(meta,o.optical_state,s.exterior_index,o.ice_material);
        meta<<",\n"
          <<"\"symmetry\":\"D6h_O3_full_24_scalar_unpolarized\",\n"
          <<"\"incident_reconstruction\":\"bilinear_theta_phi_in_spherical_incident_frame\",\n"
          <<"\"incident_formula\":\"k=(sin(t)cos(p),cos(t),sin(t)sin(p));e0=(cos(t)cos(p),-sin(t),cos(t)sin(p));e1=k_cross_e0\",\n"
          <<"\"area_quadrature\":{\"rule\":\"composite_equal_area_GL2xGL3\",\"subdivisions\":"<<s.area_subdivisions<<",\"randomness\":false},\n"
          <<"\"transport\":{\"path_rule\":\"R_and_T_R_power_m_T\",\"polarization\":\"two_input_JonesFlux\",\"interference\":false,\"diffraction\":false,\n"
          <<"\"tail_tolerance\":"<<s.tail_tolerance<<",\"max_internal_hits\":"<<s.max_internal_hits
          <<",\"max_unresolved_fraction\":"<<s.max_unresolved_fraction<<",\"max_balance_error\":"<<s.max_balance_error<<"},\n"
          <<"\"pole_policy\":\"axis_mass_uniform_in_azimuth_in_histogram_end_cap\",\n\"nodes\":[\n";
        bool first=true;std::uint64_t total_rays=0;
        for(std::uint32_t i=0;i<ni;++i)for(std::uint32_t j=0;j<nj;++j)
        {
            const auto n=atlas::make_node(s,i,j);
#if defined(ICE_ATLAS_HAS_CUDA)
            const auto r=tracer.trace(s,n);
#else
            const auto r=atlas::trace_cpu(s,n);
#endif
            mass.append(r.mass_mm2);total_rays+=r.audit.rays;
            if(first)
            {
                const std::array<std::uint64_t,1> es{std::uint64_t(n.u_edges.size())};
                rainbow::NpyFloat64Writer edges(directory.staging()/"u_edges.npy",es);edges.append(n.u_edges);edges.finish();
            }
            if(!first)meta<<",\n";
            first=false;
            meta<<"{\"i\":"<<i<<",\"j\":"<<j<<",\"theta\":"<<n.theta<<",\"phi\":"<<n.phi<<",\"ki\":";
            vec(meta,n.trace.incident_direction);meta<<",\"e0\":";vec(meta,n.grid.frame.e0);
            meta<<",\"projected_area_mm2\":"<<r.projected_area_mm2<<",\"escaped_mm2\":"<<r.escaped_mm2
              <<",\"unresolved_mm2\":"<<r.unresolved_mm2<<",\"histogram_relative_error\":"<<r.histogram_relative_error
              <<",\"rays\":"<<r.audit.rays<<",\"outputs\":"<<r.audit.outputs
              <<",\"max_balance_error\":"<<r.audit.maximum_balance_error
              <<",\"max_interface_balance_error\":"<<r.audit.maximum_interface_balance_error<<"}";
            std::cout<<"atlas_node="<<i<<','<<j<<" rays="<<r.audit.rays<<" residual="<<r.unresolved_mm2/r.projected_area_mm2<<'\n';
        }
        mass.finish();meta<<"],\n\"total_rays\":"<<total_rays<<",\"numerical_convergence_certified\":false\n}\n";meta.close();
        directory.commit();std::cout<<"atlas_complete nodes="<<ni*nj<<" rays="<<total_rays<<'\n';return 0;
    }
    catch(...){std::cerr<<"Atlas generation failed. Incomplete .part retained; no ray omission or fallback.\n";throw;}
}
int safe_run(int argc,char** argv){try{return run(argc,argv);}catch(const std::exception& e){std::cerr<<"ice_precompute_atlas: "<<e.what()<<'\n';return 2;}}
}
#if defined(_WIN32)
int wmain(int argc,wchar_t** argv)
{
    std::vector<std::string> strings;for(int i=0;i<argc;++i){auto u=std::filesystem::path(argv[i]).u8string();strings.emplace_back(reinterpret_cast<const char*>(u.data()),u.size());}
    std::vector<char*> args;for(auto& s:strings)args.push_back(s.data());return safe_run(argc,args.data());
}
#else
int main(int argc,char** argv){return safe_run(argc,argv);}
#endif
