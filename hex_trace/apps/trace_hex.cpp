#include <ice_crystal/hex_trace.hpp>
#if defined(ICE_HEX_HAS_CUDA)
#include <ice_crystal/hex_trace_cuda.hpp>
#endif
#include <algorithm>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>

namespace
{
using namespace iceCrystal;
struct Options
{
    HexPrism prism{};Vec3 incident{};
    double exterior_index=1,index=0,wavelength_nm=0,tail=1e-12;
    std::uint64_t samples=0,seed=0;
    unsigned hits=256;std::size_t batch=4096,max_outputs=4*1024*1024;
    int device=0;
#if defined(ICE_HEX_HAS_CUDA)
    std::string backend="cuda";
#else
    std::string backend="cpu";
#endif
    std::filesystem::path out,module;
};
double real(std::string_view s)
{
    double v=0;auto r=std::from_chars(s.data(),s.data()+s.size(),v);
    if(r.ec!=std::errc{}||r.ptr!=s.data()+s.size()||!finite_value(v))throw std::invalid_argument("Invalid finite real: "+std::string(s));
    return v;
}
std::uint64_t integer(std::string_view s)
{
    std::uint64_t v=0;auto r=std::from_chars(s.data(),s.data()+s.size(),v);
    if(r.ec!=std::errc{}||r.ptr!=s.data()+s.size())throw std::invalid_argument("Invalid nonnegative integer: "+std::string(s));
    return v;
}
template<class T> T bounded_integer(std::string_view s)
{
    auto n=integer(s);
    if(n>static_cast<std::uint64_t>((std::numeric_limits<T>::max)()))throw std::out_of_range("Integer option exceeds its type.");
    return static_cast<T>(n);
}
void help()
{
    std::cout<<"Hex-prism M1: fixed-orientation, flux-polarized geometrical tracing\n"
        "Required: --out DIR --radius-mm R --length-mm L --ki X Y Z\n"
        "          --ior N --wavelength-nm LAMBDA --samples N\n"
        "Options:  --external-ior 1 --seed 0 --batch-size 4096\n"
        "          --max-internal-hits 256 --tail-tolerance 1e-12\n"
        "          --max-output-records 4194304 --backend cpu|cuda --device 0\n"
        "          --module PATH_TO_FATBIN\n"
        "R is circumradius, L is FULL prism length. Axial direction is +y.\n"
        "IOR is caller-supplied at the stated vacuum wavelength; NOT an ice spectrum.\n"
        "No orientation averaging, CDF, absorption, birefringence, diffraction,\n"
        "or inter-path interference is computed by this M1 executable.\n"
        "An unaccepted result stays in DIR.part and returns exit code 2.\n";
}
std::filesystem::path path_from_utf8(std::string_view bytes)
{
    std::u8string utf8;utf8.reserve(bytes.size());
    for(unsigned char c:bytes)utf8.push_back(static_cast<char8_t>(c));
    return std::filesystem::path(utf8);
}
Options parse(int argc,char** argv)
{
    Options o;
    o.module=std::filesystem::absolute(path_from_utf8(argv[0])).parent_path()/"modules"/"hex_trace.fatbin";
    std::set<std::string> seen;
    for(int i=1;i<argc;++i)
    {
        const std::string key=argv[i];
        if(!seen.insert(key).second)throw std::invalid_argument("Duplicate option: "+key);
        auto next=[&]() -> std::string_view {if(++i>=argc)throw std::invalid_argument("Missing value for "+key);return argv[i];};
        if(key=="--out")o.out=path_from_utf8(next());
        else if(key=="--module")o.module=path_from_utf8(next());
        else if(key=="--ki"){o.incident.x=real(next());o.incident.y=real(next());o.incident.z=real(next());}
        else if(key=="--radius-mm")o.prism.circumradius_mm=real(next());
        else if(key=="--length-mm")o.prism.length_mm=real(next());
        else if(key=="--ior")o.index=real(next());
        else if(key=="--external-ior")o.exterior_index=real(next());
        else if(key=="--wavelength-nm")o.wavelength_nm=real(next());
        else if(key=="--tail-tolerance")o.tail=real(next());
        else if(key=="--samples")o.samples=integer(next());
        else if(key=="--seed")o.seed=integer(next());
        else if(key=="--max-internal-hits")o.hits=bounded_integer<unsigned>(next());
        else if(key=="--batch-size")o.batch=bounded_integer<std::size_t>(next());
        else if(key=="--max-output-records")o.max_outputs=bounded_integer<std::size_t>(next());
        else if(key=="--device")o.device=bounded_integer<int>(next());
        else if(key=="--backend")o.backend=std::string(next());
        else throw std::invalid_argument("Unknown option: "+key);
    }
    for(const char* key:{"--out","--radius-mm","--length-mm","--ki","--ior","--wavelength-nm","--samples"})
        if(!seen.contains(key))throw std::invalid_argument(std::string("Required option: ")+key);
    if(!(o.wavelength_nm>0)||!o.batch||!o.max_outputs||o.out.empty())throw std::invalid_argument("Nonpositive wavelength/batch/budget or empty output.");
    if(o.backend!="cpu"&&o.backend!="cuda")throw std::invalid_argument("Backend must be cpu or cuda.");
#if !defined(ICE_HEX_HAS_CUDA)
    if(o.backend!="cpu")throw std::invalid_argument("This is a CPU reference build; CUDA is not compiled. No fallback.");
#endif
    return o;
}
std::ofstream file(const std::filesystem::path& p)
{
    std::ofstream f(p,std::ios::binary|std::ios::trunc);
    if(!f)throw std::runtime_error("Cannot open output: "+p.string());
    f.exceptions(std::ios::badbit|std::ios::failbit);f.imbue(std::locale::classic());f<<std::setprecision(17);return f;
}
void xyz(std::ostream& f,Vec3 v){f<<v.x<<','<<v.y<<','<<v.z;}
void metadata(const std::filesystem::path& p,const Options& o,const TraceSettings& c,const TraceSummary& s)
{
    auto f=file(p);const double area=c.prism.projected_area_mm2(c.incident_direction);
    const bool complete=s.accepted&&s.rays==o.samples;
    f<<"{\n  \"schema\": \"ice.hex_trace.points.v1\",\n  \"generator_version\": \"ice.hex_trace.m1.v1\",\n"
     <<"  \"complete\": "<<(complete?"true":"false")<<",\n  \"quality_accepted\": "<<(s.accepted?"true":"false")<<",\n"
     <<"  \"phase_cdf_record\": false,\n  \"backend\": \""<<o.backend<<"\",\n"
     <<"  \"model\": \"nonabsorbing_isotropic_geometrical_optics_incoherent_paths\",\n"
     <<"  \"polarization\": \"two_input_flux_Jones_response_unpolarized_incident\",\n"
     <<"  \"tir_relative_phase_retained\": true,\n  \"inter_path_interference\": false,\n"
     <<"  \"diffraction\": false,\n  \"birefringence\": false,\n  \"absorption\": false,\n"
     <<"  \"index_model\": \"caller_supplied_real_index_at_wavelength\",\n"
     <<"  \"wavelength_vacuum_nm\": "<<o.wavelength_nm<<",\n"
     <<"  \"exterior_index\": "<<c.exterior_index<<",\n  \"interior_index\": "<<c.interior_index<<",\n"
     <<"  \"circumradius_mm\": "<<c.prism.circumradius_mm<<",\n  \"full_length_mm\": "<<c.prism.length_mm<<",\n"
     <<"  \"direction_convention\": \"physical_propagation\",\n  \"body_long_axis\": [0,1,0],\n"
     <<"  \"incident_direction\": [";xyz(f,c.incident_direction);f<<"],\n"
     <<"  \"incident_frame_e0\": [";xyz(f,c.incident_frame.e0);f<<"],\n"
     <<"  \"incident_frame_e1\": [";xyz(f,c.incident_frame.e1);f<<"],\n"
     <<"  \"projected_area_mm2\": "<<area<<",\n  \"planned_samples\": "<<o.samples<<",\n"
     <<"  \"processed_samples\": "<<s.rays<<",\n  \"seed\": "<<o.seed<<",\n"
     <<"  \"sample_generator\": \"counter_splitmix64_open52_v1\",\n"
     <<"  \"batch_size\": "<<o.batch<<",\n  \"max_output_records_per_batch\": "<<o.max_outputs<<",\n"
     <<"  \"max_internal_hits\": "<<o.hits<<",\n  \"residual_power_tolerance_per_ray\": "<<o.tail<<",\n"
     <<"  \"geometry_fp64_tolerance_epsilon_multiplier\": 128,\n"
     <<"  \"geometry_tolerance_is_rigorous_bound\": false,\n"
     <<"  \"outgoing_count\": "<<s.outputs<<",\n  \"escaped_power_sum\": "<<s.escaped_power_sum<<",\n"
     <<"  \"unresolved_power_sum\": "<<s.unresolved_power_sum<<",\n"
     <<"  \"maximum_balance_error\": "<<s.maximum_balance_error<<",\n"
     <<"  \"maximum_interface_balance_error\": "<<s.maximum_interface_balance_error<<",\n"
     <<"  \"output_weight\": \"projected_area_mm2 / planned_samples * power_fraction\",\n"
     <<"  \"renormalized_to_escaped_power\": false,\n  \"status_counts\": {";
    for(unsigned i=0;i<8;++i)f<<(i?",":"")<<"\n    \""<<status_name(static_cast<TraceStatus>(i))<<"\": "<<s.status_counts[i];
    f<<"\n  }\n}\n";f.close();
}
int run(const Options& o)
{
    auto c=make_trace_settings(o.prism,o.incident,o.exterior_index,o.index,o.samples,o.seed,o.hits,o.tail);
    auto part=o.out;part+=".part";
    if(std::filesystem::exists(o.out)||std::filesystem::exists(part))throw std::runtime_error("Output or .part already exists; refusing to overwrite.");
#if defined(ICE_HEX_HAS_CUDA)
    std::unique_ptr<rainbow::CudaContext> context;
    std::unique_ptr<HexPrismTracer> tracer;
    if(o.backend=="cuda")
    {
        context=std::make_unique<rainbow::CudaContext>(o.device);
        tracer=std::make_unique<HexPrismTracer>(*context);tracer->load_module(o.module);
    }
#endif
    std::filesystem::create_directories(part);
    auto out=file(part/"outgoing.csv"),rays=file(part/"rays.csv");
    out<<"sample_id,kind,entry_face,exit_face,internal_reflections,x_mm,y_mm,z_mm,dx,dy,dz,"
        "e0x,e0y,e0z,e1x,e1y,e1z,power_fraction,weight_mm2,internal_length_mm,"
        "j00_re,j00_im,j10_re,j10_im,j01_re,j01_im,j11_re,j11_im\n";
    rays<<"sample_id,status,output_count,internal_hits,tir_count,escaped_power,unresolved_power,balance_error,max_interface_error\n";
    TraceSummary total{};const double area_per_ray=c.prism.projected_area_mm2(c.incident_direction)/double(o.samples);
    for(std::uint64_t first=0;first<o.samples;)
    {
        const auto n=static_cast<std::size_t>((std::min)(std::uint64_t(o.batch),o.samples-first));
        TraceBatch b{};
#if defined(ICE_HEX_HAS_CUDA)
        if(tracer){tracer->trace(c,first,n,o.max_outputs);b=tracer->download();}
        else
#endif
        b=trace_batch_cpu(c,first,n,o.max_outputs);
        for(const auto& s:b.outgoing)
        {
            const double weight=area_per_ray*s.power_fraction;
            if(!finite_value(weight)||(s.power_fraction>0&&!(weight>0)))
                throw std::runtime_error("Output area weight over/underflow; refusing a completed point record.");
            out<<s.incident_sample_id<<','<<static_cast<unsigned>(s.kind)<<','<<s.entry_face<<','<<s.exit_face<<','<<s.internal_reflections<<',';
            xyz(out,s.position_mm);out<<',';xyz(out,s.direction);out<<',';xyz(out,s.frame.e0);out<<',';xyz(out,s.frame.e1);
            out<<','<<s.power_fraction<<','<<weight<<','<<s.internal_length_mm;
            for(const auto& col:s.response.column)out<<','<<col.x.re<<','<<col.x.im<<','<<col.y.re<<','<<col.y.im;
            out<<'\n';
        }
        for(std::size_t i=0;i<b.audits.size();++i)
        {
            const auto& a=b.audits[i];
            rays<<first+i<<','<<status_name(a.status)<<','<<a.output_count<<','<<a.internal_hits<<','<<a.tir_count<<','
                <<a.escaped_power<<','<<a.unresolved_power<<','<<a.balance_error<<','<<a.max_interface_balance_error<<'\n';
        }
        const auto summary=summarize(b.audits);merge_summary(total,summary);first+=n;
        if(!summary.accepted)break;
    }
    out.close();rays.close();metadata(part/"trace_metadata.json",o,c,total);
    std::cout<<std::setprecision(17)<<"backend="<<o.backend<<" rays="<<total.rays<<" outputs="<<total.outputs
        <<" escaped_mean="<<total.escaped_power_sum/double(total.rays)
        <<" unresolved_mean="<<total.unresolved_power_sum/double(total.rays)
        <<" max_balance_error="<<total.maximum_balance_error<<'\n';
    if(!total.accepted||total.rays!=o.samples)
    {std::cerr<<"Trace rejected; diagnostics remain in "<<part<<". No normalization or success commit.\n";return 2;}
    std::filesystem::rename(part,o.out);return 0;
}
}
int run_main(int argc,char** argv)
{
    try
    {
        if(argc==2 && std::string_view(argv[1])=="--help"){help();return 0;}
        return run(parse(argc,argv));
    }
    catch(const std::exception& e){std::cerr<<"ice_trace_hex: "<<e.what()<<'\n';return 1;}
}

#if defined(_WIN32)
int wmain(int argc,wchar_t** argv)
{
    try
    {
        std::vector<std::string> storage;storage.reserve(static_cast<std::size_t>(argc));
        for(int i=0;i<argc;++i)
        {
            const auto u=std::filesystem::path(argv[i]).u8string();
            storage.emplace_back(reinterpret_cast<const char*>(u.data()),u.size());
        }
        std::vector<char*> args;args.reserve(storage.size());
        for(auto& a:storage)args.push_back(a.data());
        return run_main(argc,args.data());
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
#else
int main(int argc,char** argv){return run_main(argc,argv);}
#endif
