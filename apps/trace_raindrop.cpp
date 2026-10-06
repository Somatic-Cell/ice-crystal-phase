#include <rainbow/raindrop_tracer.hpp>
#include <rainbow/patch_accel.hpp>
#include <rainbow/patch_query.hpp>
#include <rainbow/patch_optics.hpp>
#include <rainbow/patch_failure_report.hpp>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <span>
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>

namespace
{
struct TraceCommandLine
{
    std::filesystem::path module_path,output_path;
    std::filesystem::path patch_module_path,patch_csv_path;
    std::filesystem::path query_module_path,query_csv_path,query_hits_csv_path;
    std::filesystem::path optics_module_path,optics_csv_path;
    std::filesystem::path wave_csv_path;
    std::filesystem::path patch_failure_report_path;
    rainbow::WaveOpticsSettings wave_settings{};
    bool focal_offsets_given=false, diffraction_option_given=false;
    std::uint32_t query_theta=90,query_phi=180;
    rainbow::RaindropSettings settings;

    template<class Char> static TraceCommandLine parse(const int argc,Char* argv[])
    {
        if(argc<3) throw std::invalid_argument(
            "Usage: rainbow_trace <raindrop_trace.optixir> <vertices.csv> [--radius-mm value] "
            "[--grid count] [--inclination-deg value] [--azimuth-deg value] "
            "[--wavelength-nm value --ior value] [--polarization x|y] [--sphere] "
            "[--patch-module patch_build.fatbin] [--patch-csv patches.csv] "
            "[--query-module patch_query.optixir --query-csv queries.csv] "
            "[--query-hits-csv hits.csv] [--query-theta N --query-phi N] "
            "[--optics-module patch_optics.fatbin --optics-csv optics.csv] "
            "[--wave-csv wave.csv --focal-offsets R,TT,TRT,TRRT] "
            "[--diffraction-sigma-deg value] [--diffraction-contrast value] "
            "[--patch-failure-report NEW_report.json]");
        TraceCommandLine command;
        command.module_path=std::filesystem::path(argv[1]);
        command.output_path=std::filesystem::path(argv[2]);
        float inclination=0.0f,azimuth=0.0f;
        bool wavelength_given=false,index_given=false;
        const auto ascii=[](const Char* input)
        {
            std::string result;
            for(;*input;++input)
            {
                const unsigned long c=static_cast<unsigned long>(*input);
                if(c>127)throw std::invalid_argument("Option names must be ASCII.");
                result.push_back(static_cast<char>(c));
            }
            return result;
        };
        for(int i=3;i<argc;++i)
        {
            const std::string option=ascii(argv[i]);
            if(option=="--sphere"){command.settings.force_sphere=true;continue;}
            if(++i>=argc)throw std::invalid_argument("Missing option value.");
            // Path values keep the native character type (including Japanese paths).
            if(option=="--patch-module")
            {command.patch_module_path=std::filesystem::path(argv[i]);continue;}
            if(option=="--patch-csv")
            {command.patch_csv_path=std::filesystem::path(argv[i]);continue;}
            if(option=="--query-module")
            {command.query_module_path=std::filesystem::path(argv[i]);continue;}
            if(option=="--query-csv")
            {command.query_csv_path=std::filesystem::path(argv[i]);continue;}
            if(option=="--query-hits-csv")
            {command.query_hits_csv_path=std::filesystem::path(argv[i]);continue;}
            if(option=="--optics-module")
            {command.optics_module_path=std::filesystem::path(argv[i]);continue;}
            if(option=="--optics-csv")
            {command.optics_csv_path=std::filesystem::path(argv[i]);continue;}
            if(option=="--patch-failure-report")
            {command.patch_failure_report_path=std::filesystem::path(argv[i]);continue;}
            if(option=="--wave-csv")
            {command.wave_csv_path=std::filesystem::path(argv[i]);continue;}
            const std::string value=ascii(argv[i]);
            if(option=="--focal-offsets")
            {
                // Require an explicit baseline convention. Not inferred from a
                // signed patch area or silently calibrated against Mie curves.
                if(value.size()!=7 || value[1]!=',' || value[3]!=',' || value[5]!=',')
                    throw std::invalid_argument("focal-offsets must be four digits in 0..3, e.g. 0,0,0,0.");
                for(unsigned j=0;j<4;++j)
                {
                    const char c=value[2*j];
                    if(c<'0'||c>'3')throw std::invalid_argument("Each focal offset is 0..3 quarter-turns.");
                    command.wave_settings.focal_quarter_turn_offsets[j]=static_cast<std::uint32_t>(c-'0');
                }
                command.focal_offsets_given=true;continue;
            }
            if(option=="--diffraction-sigma-deg" || option=="--diffraction-contrast")
            {
                std::size_t used=0;const double x=std::stod(value,&used);
                if(used!=value.size() || !std::isfinite(x))throw std::invalid_argument("Expected a finite diffraction option.");
                if(option=="--diffraction-sigma-deg")
                {
                    if(!(x>0))throw std::invalid_argument("Explicit diffraction sigma must be positive.");
                    command.wave_settings.primary_sigma_degrees=x;
                }
                else
                {
                    if(!(x>=0 && x<1))throw std::invalid_argument("Diffraction contrast must be in [0,1).");
                    command.wave_settings.transition_contrast=x;
                }
                command.diffraction_option_given=true;continue;
            }
            if(option=="--polarization")
            {
                if(value=="x")command.settings.incident_field={{1,0},{0,0}};
                else if(value=="y")command.settings.incident_field={{0,0},{1,0}};
                else throw std::invalid_argument("polarization must be x or y (separate coherent input states).");
                continue;
            }
            if(option=="--grid" || option=="--query-theta" || option=="--query-phi")
            {
                if(value.empty()||value.find_first_not_of("0123456789")!=std::string::npos)
                    throw std::invalid_argument("grid must be an unsigned integer.");
                std::size_t used=0;const unsigned long count=std::stoul(value,&used);
                if(used!=value.size()||count>32767)throw std::invalid_argument("grid is out of range.");
                if(count==0)throw std::invalid_argument("Grid count must be positive.");
                if(option=="--query-theta")command.query_theta=static_cast<std::uint32_t>(count);
                else if(option=="--query-phi")command.query_phi=static_cast<std::uint32_t>(count);
                else command.settings.grid_width=command.settings.grid_height=static_cast<std::uint32_t>(count);
                continue;
            }
            std::size_t used=0;const float number=std::stof(value,&used);
            if(used!=value.size()||!std::isfinite(number))throw std::invalid_argument("Expected a finite numerical option.");
            if(option=="--radius-mm")command.settings.radius_mm=number;
            else if(option=="--wavelength-nm"){command.settings.wavelength_nm=number;wavelength_given=true;}
            else if(option=="--ior"){command.settings.interior_index=number;index_given=true;}
            else if(option=="--inclination-deg")inclination=number;
            else if(option=="--azimuth-deg")azimuth=number;
            else throw std::invalid_argument("Unknown option: "+option);
        }
        if(!command.patch_csv_path.empty() && command.patch_module_path.empty())
            throw std::invalid_argument("--patch-csv requires --patch-module.");
        if(!command.patch_csv_path.empty()
           && std::filesystem::absolute(command.patch_csv_path).lexically_normal()
              ==std::filesystem::absolute(command.output_path).lexically_normal())
            throw std::invalid_argument("Vertex CSV and patch CSV must have different paths.");
        const bool wants_query=!command.query_module_path.empty();
        if(wants_query && (command.patch_module_path.empty() || command.query_csv_path.empty()))
            throw std::invalid_argument("--query-module requires --patch-module and --query-csv.");
        if(!wants_query && (!command.query_csv_path.empty() || !command.query_hits_csv_path.empty()))
            throw std::invalid_argument("Query output requires --query-module.");
        const bool wants_optics=!command.optics_module_path.empty();
        if(wants_optics != !command.optics_csv_path.empty())
            throw std::invalid_argument("--optics-module and --optics-csv must be specified together.");
        if(wants_optics && !wants_query)
            throw std::invalid_argument("Optical evaluation requires --query-module and --query-csv.");
        const bool wants_wave=!command.wave_csv_path.empty();
        if(wants_wave && (!wants_optics || !command.focal_offsets_given))
            throw std::invalid_argument("--wave-csv requires optical evaluation and explicit --focal-offsets R,TT,TRT,TRRT.");
        if(!wants_wave && (command.focal_offsets_given || command.diffraction_option_given))
            throw std::invalid_argument("Focal/diffraction options require --wave-csv.");
        if(!command.patch_failure_report_path.empty() && !wants_wave)
            throw std::invalid_argument("--patch-failure-report requires --wave-csv.");
        if(!command.patch_failure_report_path.empty())
        {
            auto partial=command.patch_failure_report_path;partial+=".part";
            if(std::filesystem::exists(command.patch_failure_report_path)||std::filesystem::exists(partial))
                throw std::invalid_argument("Choose a new --patch-failure-report path (report or .part already exists).");
        }
        const std::array<std::filesystem::path,7> outputs={command.output_path,command.patch_csv_path,
            command.query_csv_path,command.query_hits_csv_path,command.optics_csv_path,command.wave_csv_path,command.patch_failure_report_path};
        for(std::size_t a=0;a<outputs.size();++a)for(std::size_t b=a+1;b<outputs.size();++b)
        {
            if(outputs[a].empty() || outputs[b].empty())continue;
            const auto pa=std::filesystem::absolute(outputs[a]).lexically_normal();
            const auto pb=std::filesystem::absolute(outputs[b]).lexically_normal();
            if(pa==pb || (std::filesystem::exists(pa) && std::filesystem::exists(pb)
                         && std::filesystem::equivalent(pa,pb)))
                throw std::invalid_argument("Output paths must be distinct.");
        }
        if(wavelength_given&&!index_given)
            throw std::invalid_argument("Specify --ior together with --wavelength-nm; no water dispersion fit is silently assumed.");
        if(std::abs(inclination)>90||std::abs(azimuth)>360)
            throw std::invalid_argument("inclination must be in [-90,90], azimuth in [-360,360] degrees.");
        constexpr float radians=0.01745329251994329577f;
        const float e=inclination*radians,a=azimuth*radians;
        // +y は鉛直上向き．inclination>0 は下降する入射光，azimuth=0 は +x 方向．
        command.settings.incident_direction={std::cos(e)*std::cos(a),-std::sin(e),std::cos(e)*std::sin(a)};
        static_cast<void>(command.settings.make_config());
        return command;
    }
};

template<class Char> int run(const int argc,Char* argv[])
{
    try
    {
        const auto command = TraceCommandLine::parse(argc, argv);

        rainbow::CudaContext cuda_context{0};
        rainbow::RaindropTracer tracer(cuda_context);

        tracer.create_pipeline(command.module_path);
        tracer.trace(command.settings);
        // 診断があっても CSV に保存する．数値失敗は成功とせず exit code 1 で通知する．
        const auto vertices=tracer.download_vertices();
        tracer.write_csv(command.output_path,std::span<const rainbow::OutgoingVertex>{vertices});
        const auto count=tracer.config().vertex_count();
        constexpr const char* names[]={"R","TT","TRT","TRRT"};
        std::size_t errors=0;
        for(unsigned family=0;family<4;++family)
        {
            std::array<std::size_t,7> histogram{};
            std::size_t fallback=0;
            for(std::uint32_t i=0;i<count;++i)
            {
                const auto& v=vertices[std::size_t(family)*count+i];
                const auto status=static_cast<unsigned>(v.status);
                if(status>=histogram.size())throw std::runtime_error("Corrupt vertex status.");
                ++histogram[status];fallback+=bool(v.diagnostics&rainbow::IntersectionFallback);
            }
            const auto failed=histogram[0]+histogram[4]+histogram[5]+histogram[6];errors+=failed;
            std::cout<<names[family]<<": valid="<<histogram[1]<<", miss="<<histogram[2]
                <<", TIR="<<histogram[3]<<", numerical_failures="<<failed<<", fallback="<<fallback<<'\n';
        }
        std::cout<<"Saved "<<command.output_path<<" (outgoing vertices, NOT a phase-function LUT).\n";
        if(errors){std::cerr<<"Unresolved vertices remain; inspect status/diagnostics before patch construction.\n";return EXIT_FAILURE;}

        if(!command.patch_module_path.empty())
        {
            // The tracer and its vertices outlive PatchAccel. No re-upload of
            // the CSV/readback vertices: build() consumes the existing GPU buffer.
            rainbow::PatchAccel patches(cuda_context,tracer.optix_context());
            patches.load_module(command.patch_module_path);
            patches.build(tracer.vertices(),tracer.config());
            const auto& stats=patches.statistics();
            using Status=rainbow::PatchCellStatus;
            std::cout<<"Patch geometry: logical_cells="<<stats.logical_cell_count
                <<", stored_patches="<<stats.patch_count
                <<", regular_positive="<<stats.count(Status::RegularPositive)
                <<", regular_negative="<<stats.count(Status::RegularNegative)
                <<", needs_refinement="<<stats.count(Status::NeedsRefinement)
                <<", missing_corners="<<stats.count(Status::MissingCorners)
                <<", no_outgoing_corners="<<stats.count(Status::NoOutgoingCorners)
                <<", invalid_vertices="<<stats.count(Status::InvalidVertex)
                <<", invalid_geometry="<<stats.count(Status::InvalidGeometry)<<'\n';
            std::cout<<"Patch GAS: bytes="<<stats.gas_byte_size
                <<", nonempty="<<(patches.handle()!=0)<<'\n';
            std::cout<<"This is geometric assembly, NOT optical completion: "
                <<"boundary cells and nonregular maps are retained for later resolution.\n";
            if(!command.patch_csv_path.empty())
            {
                patches.write_csv(command.patch_csv_path);
                std::cout<<"Saved "<<command.patch_csv_path<<" (all logical cell statuses).\n";
            }
            if(!command.query_module_path.empty())
            {
                // patches と tracer は query より長く生存する．GAS は再構築しない．
                rainbow::PatchQuery query(cuda_context,tracer.optix_context());
                query.create_pipeline(command.query_module_path);
                query.query_grid(patches,tracer.config(),command.query_theta,command.query_phi);
                query.write_csv(command.query_csv_path);
                if(!command.query_hits_csv_path.empty())query.write_hits_csv(command.query_hits_csv_path);
                const auto& q=query.statistics();
                std::cout<<"Patch query: directions="<<q.directions<<", nonempty="<<q.nonempty_directions
                    <<", intersections="<<q.hits<<", max_hits="<<q.max_hits
                    <<", unresolved_or_error_directions="<<q.error_directions
                    <<", refinement_directions="<<q.refinement_directions
                    <<", boundary_directions="<<q.boundary_directions
                    <<", fp64_directions="<<q.fp64_directions<<'\n';
                std::cout<<"Saved "<<command.query_csv_path<<" (geometric hit counts, NOT a phase function).\n";
                std::cout<<"Missing source cells / boundary ownership / nonregular optical evaluation remain pending.\n";
                bool failed=q.error_directions!=0;
                if(!command.optics_module_path.empty())
                {
                    // Consume the existing GPU arrays before query/patches/tracer
                    // are closed. No new OptiX pipeline and no CSV round trip.
                    rainbow::PatchOptics optics(cuda_context);
                    optics.load_module(command.optics_module_path);
                    if(command.wave_csv_path.empty()) optics.evaluate(patches,query,tracer.config());
                    else optics.evaluate_wave(patches,query,tracer.config(),command.wave_settings);
                    // Read-only capture BEFORE final CSV I/O and before owners close.
                    // The report does not alter any numerical output or error mask.
                    if(!command.patch_failure_report_path.empty())
                    {
                        const auto report=rainbow::PatchFailureReport::capture(
                            cuda_context,patches,query,optics,tracer.config(),command.wave_settings);
                        report.write_json(command.patch_failure_report_path);
                        std::cout<<"Patch failure report: directions="<<report.snapshot().directions.size()
                            <<", patches="<<report.snapshot().patches.size()
                            <<", sparse_read_bytes="<<report.snapshot().gpu_read_bytes<<'\n';
                        std::cout<<"Saved "<<command.patch_failure_report_path<<" (read-only witnesses, NOT repaired values).\n";
                    }
                    optics.write_csv(command.optics_csv_path);
                    bool wave_failed=false;
                    if(!command.wave_csv_path.empty())
                    {
                        optics.write_wave_csv(command.wave_csv_path);
                        const auto& w=optics.wave_statistics();
                        std::cout<<"Wave optics: focal_errors="<<w.focal_errors<<", focal_pending="<<w.focal_pending
                            <<", corrected_hits="<<w.corrected_hits<<", extra_quarter_turn_hits="<<w.extra_quarter_turn_hits
                            <<", primary_transitions="<<w.primary_transitions<<", secondary_transitions="<<w.secondary_transitions
                            <<", unknown_transitions="<<w.unknown_transitions<<", filtered_directions="<<w.filtered_directions
                            <<", unavailable_directions="<<w.unavailable_directions<<", underresolved_directions="<<w.underresolved_directions<<'\n';
                        std::cout<<"Saved "<<command.wave_csv_path<<" (focal and local Gaussian diffraction; source coverage remains uncertified).\n";
                        if(w.underresolved_directions || w.filtered_directions==0u)
                            std::cerr<<"Warning: verify angular resolution relative to the tabulated sigma; no bandwidth enlargement is performed.\n";
                        wave_failed=w.focal_errors!=0u;
                        for(const auto& d:optics.host_diffraction_results())
                            wave_failed=wave_failed||((d.flags&rainbow::DiffractionInvalidData)!=0u);
                    }
                    const auto& optical=optics.statistics();
                    std::cout<<"Patch optics: directions="<<optical.directions
                        <<", known_hits_complete_directions="<<optical.known_hits_complete_directions
                        <<", pending_directions="<<optical.pending_directions
                        <<", error_directions="<<optical.error_directions
                        <<", evaluated_hits="<<optical.evaluated_hits
                        <<", rejected_hits="<<optical.rejected_hits<<'\n';
                    std::cout<<"Saved "<<command.optics_csv_path
                        <<" (regular partial optical densities, NOT a complete phase function).\n";
                    if(command.wave_csv_path.empty())
                        std::cout<<"Focal-line phase and diffraction are not enabled in this run.\n";
                    std::cout<<"Source coverage and boundary ownership remain uncertified.\n";
                    failed=failed || optical.error_directions!=0 || wave_failed;
                    optics.close();
                }
                query.close();
                if(failed)throw std::runtime_error("Query/optical numerical errors remain; see CSV flags.");
            }
            patches.close();
        }
        tracer.close();
        return EXIT_SUCCESS;
    }
    catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<'\n';return EXIT_FAILURE;}
}
}
#if defined(_WIN32)
int wmain(const int argc,wchar_t* argv[]){return run(argc,argv);}
#else
int main(const int argc,char* argv[]){return run(argc,argv);}
#endif
