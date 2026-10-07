#include <rainbow/raindrop_tracer.hpp>
#include <rainbow/patch_accel.hpp>
#include <rainbow/patch_query.hpp>
#include <rainbow/patch_optics.hpp>
#include <rainbow/phase_cdf.hpp>
#include <rainbow/phase_cdf_math.hpp>
#include <rainbow/npy_writer.hpp>
#include <rainbow/folded_patch_data.hpp>

#include <array>
#include <charconv>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <stdexcept>
#include <string>
#include <string_view>

#ifndef RAINBOW_SOURCE_COMMIT
#define RAINBOW_SOURCE_COMMIT "unknown"
#endif
#ifndef RAINBOW_SOURCE_DIRTY
#define RAINBOW_SOURCE_DIRTY 1
#endif

namespace
{
using namespace rainbow;
struct Command
{
    RaindropSettings settings{};
    WaveOpticsSettings wave{};
    PhaseCdfPolicy policy{};
    std::filesystem::path modules, output;
    std::uint32_t theta=360,phi=720;
    PhaseDensityStage stage=PhaseDensityStage::Diffraction;
    int device=0;
};
double number(std::string_view s)
{
    double x=0;
    const auto r=std::from_chars(s.data(),s.data()+s.size(),x);
    if(r.ec!=std::errc{} || r.ptr!=s.data()+s.size() || !std::isfinite(x))
        throw std::invalid_argument("Invalid numeric option.");
    return x;
}
std::uint32_t positive_integer(std::string_view s)
{
    std::uint32_t x=0;
    const auto r=std::from_chars(s.data(),s.data()+s.size(),x);
    if(r.ec!=std::errc{} || r.ptr!=s.data()+s.size() || x==0)
        throw std::invalid_argument("Invalid positive integer option.");
    return x;
}
const char* stage_name(PhaseDensityStage s)
{
    switch(s)
    {
        case PhaseDensityStage::Incoherent: return "incoherent";
        case PhaseDensityStage::Path: return "path";
        case PhaseDensityStage::Focal: return "focal";
        case PhaseDensityStage::Diffraction: return "diffraction";
        default: throw std::invalid_argument("Invalid optical stage.");
    }
}
template<class Char> Command parse(int argc,Char** argv)
{
    Command c;c.settings.radius_mm=1.0f;c.settings.grid_width=c.settings.grid_height=513;
    c.modules=std::filesystem::absolute(std::filesystem::path(argv[0])).parent_path()/"modules";
    double inclination=20;
    for(int i=1;i<argc;++i)
    {
        const auto option=std::filesystem::path(argv[i]).string();
        if(option=="--sphere") {c.settings.force_sphere=true;continue;}
        if(option=="--allow-underresolved") {c.policy.allow_underresolved=true;continue;}
        if(option=="--help")
        {
            std::cout << "rainbow_generate --out DIRECTORY [--modules DIRECTORY] [--sphere]\n"
                         "  --radius-mm R --wavelength-nm L --ior N --grid VERTICES_PER_AXIS\n"
                         "  --inclination-deg A --query-theta T --query-phi P\n"
                         "  --stage diffraction|focal|incoherent|path\n"
                         "  --diffraction-sigma-deg S --focal-offsets 0,0,0,0\n"
                         "  --allow-underresolved  (explicitly accept flagged angular underresolution)\n"
                         "  --cdf-max-l1 E --cdf-max-lost-mass E --device ORDINAL\n";
            return {};
        }
        if(++i>=argc) throw std::invalid_argument("Missing option value.");
        if(option=="--out") {c.output=std::filesystem::path(argv[i]);continue;}
        if(option=="--modules") {c.modules=std::filesystem::path(argv[i]);continue;}
        const auto value=std::filesystem::path(argv[i]).string();
        if(option=="--radius-mm") c.settings.radius_mm=static_cast<float>(number(value));
        else if(option=="--wavelength-nm") c.settings.wavelength_nm=static_cast<float>(number(value));
        else if(option=="--ior") c.settings.interior_index=static_cast<float>(number(value));
        else if(option=="--grid") c.settings.grid_width=c.settings.grid_height=positive_integer(value);
        else if(option=="--query-theta") c.theta=positive_integer(value);
        else if(option=="--query-phi") c.phi=positive_integer(value);
        else if(option=="--inclination-deg") inclination=number(value);
        else if(option=="--diffraction-sigma-deg") c.wave.primary_sigma_degrees=number(value);
        else if(option=="--transition-contrast") c.wave.transition_contrast=number(value);
        else if(option=="--cdf-max-l1") c.policy.maximum_l1_error=number(value);
        else if(option=="--cdf-max-lost-mass") c.policy.maximum_lost_mass=number(value);
        else if(option=="--device")
        {
            const auto r=std::from_chars(value.data(),value.data()+value.size(),c.device);
            if(r.ec!=std::errc{} || r.ptr!=value.data()+value.size() || c.device<0)
                throw std::invalid_argument("Invalid device ordinal.");
        }
        else if(option=="--stage")
        {
            if(value=="incoherent")c.stage=PhaseDensityStage::Incoherent;
            else if(value=="path")c.stage=PhaseDensityStage::Path;
            else if(value=="focal")c.stage=PhaseDensityStage::Focal;
            else if(value=="diffraction")c.stage=PhaseDensityStage::Diffraction;
            else throw std::invalid_argument("Unknown stage.");
        }
        else if(option=="--focal-offsets")
        {
            std::size_t begin=0;
            for(unsigned k=0;k<4;++k)
            {
                const auto end=value.find(',',begin);
                if((k<3 && end==std::string::npos) || (k==3 && end!=std::string::npos))
                    throw std::invalid_argument("Expected four focal offsets.");
                const auto length=end==std::string::npos?value.size()-begin:end-begin;
                const auto part=std::string_view(value).substr(begin,length);
                auto& x=c.wave.focal_quarter_turn_offsets[k];
                const auto r=std::from_chars(part.data(),part.data()+part.size(),x);
                if(r.ec!=std::errc{} || r.ptr!=part.data()+part.size() || x>3)
                    throw std::invalid_argument("Focal offsets must be 0..3.");
                begin=end==std::string::npos?value.size():end+1;
            }
        }
        else throw std::invalid_argument("Unknown option: "+option);
    }
    if(c.output.empty()) throw std::invalid_argument("--out is required.");
    const double angle=inclination*phase_cdf_pi/180.0;
    c.settings.incident_direction={static_cast<float>(std::cos(angle)),static_cast<float>(-std::sin(angle)),0};
    static_cast<void>(c.settings.make_config());
    if(std::filesystem::exists(c.output)) throw std::runtime_error("Output already exists.");
    auto partial=c.output;partial+=".part";
    if(std::filesystem::exists(partial)) throw std::runtime_error("Previous .part directory exists; inspect it before recomputing.");
    return c;
}
using D3=std::array<double,3>;
double dot(D3 a,D3 b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
D3 normalized(D3 a)
{
    const double n=std::sqrt(dot(a,a));
    if(!(n>0 && std::isfinite(n)))throw std::runtime_error("Invalid frame.");
    for(auto& v:a)v/=n;
    return a;
}
void write_metadata(const std::filesystem::path& path,const Command& c,
    const RaindropTraceConfig& config,const PatchAccel& patches,const PhaseCdf& cdf)
{
    std::ofstream out(path,std::ios::binary|std::ios::trunc);
    if(!out)throw std::runtime_error("Cannot open metadata.json.");
    out.imbue(std::locale::classic());out<<std::setprecision(17);
    const auto& s=cdf.input_statistics();const auto& a=cdf.audit_statistics();
    auto k=normalized({config.incident_direction.x,config.incident_direction.y,config.incident_direction.z});
    D3 e{config.incident_basis_x.x,config.incident_basis_x.y,config.incident_basis_x.z};
    const double d=dot(e,k);for(unsigned i=0;i<3;++i)e[i]-=d*k[i];e=normalized(e);
    const D3 f{k[1]*e[2]-k[2]*e[1],k[2]*e[0]-k[0]*e[2],k[0]*e[1]-k[1]*e[0]};
    const FoldedPatchConfig folded{};
    out << "{\n  \"schema\": \"rainbow.phase_cdf.numpy.v1\",\n  \"complete\": true,\n"
        << "  \"dtype\": \"<f8\",\n  \"order\": \"C\",\n"
        << "  \"theta_count\": "<<c.theta<<",\n  \"phi_count\": "<<c.phi<<",\n"
        << "  \"density_measure\": \"solid_angle_sr\",\n"
        << "  \"cell_model\": \"constant_density_per_spherical_cell\",\n"
        << "  \"coordinates\": \"u=(1-cos(theta))/2; v=(phi+pi)/(2*pi)\",\n"
        << "  \"theta_definition\": \"angle_between_incident_and_outgoing_propagation_directions\",\n"
        << "  \"cdf_axis_order\": [\"phi_cell\", \"theta_edge\"],\n"
        << "  \"phi_range_rad\": ["<<-phase_cdf_pi<<','<<phase_cdf_pi<<"],\n"
        << "  \"normalization\": \"integral_p_domega_equals_one\",\n"
        << "  \"input_polarization\": \"unpolarized\",\n"
        << "  \"incident_total_intensity\": 1,\n  \"single_scattering_albedo\": 1,\n"
        << "  \"collision_model\": \"geometric_projected_area\",\n"
        << "  \"projected_area_storage\": \"separate_geometry_table_not_implemented_here\",\n"
        << "  \"stage\": \""<<stage_name(c.stage)<<"\",\n"
        << "  \"radius_mm\": "<<config.radius_mm<<",\n"
        << "  \"wavelength_nm\": "<<config.wavelength_nm<<",\n"
        << "  \"exterior_index\": "<<config.exterior_index<<",\n"
        << "  \"interior_index\": "<<config.interior_index<<",\n"
        << "  \"force_sphere\": "<<(c.settings.force_sphere?"true":"false")<<",\n"
        << "  \"incident_grid\": ["<<config.grid_width<<','<<config.grid_height<<"],\n"
        << "  \"incident_grid_half_extent\": "<<config.grid_half_extent<<",\n"
        << "  \"reference_distances\": ["<<config.reference_distance<<','<<config.outgoing_reference_distance<<"],\n"
        << "  \"query_frame_axis\": ["<<config.incident_direction.x<<','<<config.incident_direction.y<<','<<config.incident_direction.z<<"],\n"
        << "  \"query_frame_e0\": ["<<config.incident_basis_x.x<<','<<config.incident_basis_x.y<<','<<config.incident_basis_x.z<<"],\n"
        << "  \"sampling_frame_columns\": [";
    for(unsigned r=0;r<3;++r)out<<(r?",":"")<<'['<<e[r]<<','<<f[r]<<','<<k[r]<<']';
    out<<"],\n  \"shape_coefficients\": [";
    for(unsigned i=0;i<8;++i)out<<(i?",":"")<<config.shape.coefficients[i];
    out<<"],\n  \"focal_quarter_turn_offsets\": [";
    for(unsigned i=0;i<4;++i)out<<(i?",":"")<<c.wave.focal_quarter_turn_offsets[i];
    out<<"],\n  \"diffraction_primary_sigma_degrees_setting\": "<<c.wave.primary_sigma_degrees
       <<",\n  \"diffraction_zero_setting_selects_Table_II\": true,\n"
       <<"  \"diffraction_transition_contrast\": "<<c.wave.transition_contrast<<",\n"
       <<"  \"folded_branch_model\": \"finite_branch_area_ratio_v1\",\n"
       <<"  \"folded_relative_tolerance\": "<<folded.relative_tolerance<<",\n"
       <<"  \"source_commit\": \""<<RAINBOW_SOURCE_COMMIT<<"\",\n"
       <<"  \"source_dirty_at_configure\": "<<(RAINBOW_SOURCE_DIRTY?"true":"false")<<",\n"
       <<"  \"generator_version\": \"numpy_cdf_v1\",\n"
       <<"  \"quality\": {\n    \"source_coverage_certified\": false,\n"
       <<"    \"angular_convergence_certified\": false,\n"
       <<"    \"missing_source_cells\": "<<patches.statistics().count(PatchCellStatus::MissingCorners)<<",\n"
       <<"    \"invalid_values\": "<<s.invalid_values<<",\n"
       <<"    \"incomplete_directions\": "<<s.incomplete<<",\n"
       <<"    \"underresolved_directions\": "<<s.underresolved<<",\n"
       <<"    \"allow_underresolved\": "<<(c.policy.allow_underresolved?"true":"false")<<",\n"
       <<"    \"cdf_l1_mass_error\": "<<a.l1_error<<",\n"
       <<"    \"cdf_maximum_cell_mass_error\": "<<a.maximum_cell_error<<",\n"
       <<"    \"cdf_lost_positive_cells\": "<<a.lost_cells<<",\n"
       <<"    \"cdf_lost_probability_mass\": "<<a.lost_mass<<",\n"
       <<"    \"maximum_accepted_l1_mass_error\": "<<c.policy.maximum_l1_error<<",\n"
       <<"    \"maximum_accepted_lost_probability_mass\": "<<c.policy.maximum_lost_mass<<"\n  }\n}\n";
    out.flush();if(!out)throw std::runtime_error("Metadata write failed.");
    out.close();if(out.fail())throw std::runtime_error("Metadata close failed.");
}
template<class Char> int run(int argc,Char** argv)
{
    const auto c=parse(argc,argv);
    if(c.output.empty())return 0;
    const auto start=std::chrono::steady_clock::now();
    CudaContext context(c.device);
    RaindropTracer tracer(context);
    PatchAccel patches(context,tracer.optix_context());
    PatchQuery query(context,tracer.optix_context());
    PatchOptics optics(context);
    PhaseCdf cdf(context);
    tracer.create_pipeline(c.modules/"raindrop_trace.optixir");
    query.create_pipeline(c.modules/"patch_query.optixir");
    patches.load_module(c.modules/"patch_build.fatbin");
    optics.load_module(c.modules/"patch_optics.fatbin");
    cdf.load_module(c.modules/"phase_cdf.fatbin");
    tracer.trace_unpolarized(c.settings);
    cdf.validate_trace(tracer); // Small GPU reduction, no vertex CSV or bulk readback.
    patches.build(tracer);
    query.query_grid(patches,tracer.config(),c.theta,c.phi);
    if(query.statistics().error_directions)throw std::runtime_error("Patch query has numerical failures.");
    if(c.stage==PhaseDensityStage::Incoherent || c.stage==PhaseDensityStage::Path)
        optics.evaluate_device(patches,query,tracer.config());
    else optics.evaluate_wave_device(patches,query,tracer.config(),c.wave);
    cdf.build(optics,c.theta,c.phi,c.stage,c.policy);
    DatasetDirectory destination(c.output);
    cdf.write_arrays(destination.staging());
    write_metadata(destination.staging()/"metadata.json",c,tracer.config(),patches,cdf);
    destination.commit();
    const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"Saved NumPy CDF dataset: "<<c.output<<"; elapsed_seconds="<<elapsed
        <<"; CDF_L1_mass_error="<<cdf.audit_statistics().l1_error
        <<"; lost_probability_mass="<<cdf.audit_statistics().lost_mass<<'\n';
    return 0;
}
}
#if defined(_WIN32)
int wmain(int argc,wchar_t** argv)
#else
int main(int argc,char** argv)
#endif
{
    try{return run(argc,argv);}
    catch(const std::exception& e){std::cerr<<"rainbow_generate: "<<e.what()<<'\n';return 1;}
}
