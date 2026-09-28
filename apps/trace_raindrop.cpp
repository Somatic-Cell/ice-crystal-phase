#include <rainbow/raindrop_tracer.hpp>
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
    rainbow::RaindropSettings settings;

    template<class Char> static TraceCommandLine parse(const int argc,Char* argv[])
    {
        if(argc<3) throw std::invalid_argument(
            "Usage: rainbow_trace <raindrop_trace.optixir> <vertices.csv> [--radius-mm value] "
            "[--grid count] [--inclination-deg value] [--azimuth-deg value] "
            "[--wavelength-nm value --ior value] [--polarization x|y] [--sphere]");
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
            const std::string value=ascii(argv[i]);
            if(option=="--polarization")
            {
                if(value=="x")command.settings.incident_field={{1,0},{0,0}};
                else if(value=="y")command.settings.incident_field={{0,0},{1,0}};
                else throw std::invalid_argument("polarization must be x or y (separate coherent input states).");
                continue;
            }
            if(option=="--grid")
            {
                if(value.empty()||value.find_first_not_of("0123456789")!=std::string::npos)
                    throw std::invalid_argument("grid must be an unsigned integer.");
                std::size_t used=0;const unsigned long count=std::stoul(value,&used);
                if(used!=value.size()||count>32767)throw std::invalid_argument("grid is out of range.");
                command.settings.grid_width=command.settings.grid_height=static_cast<std::uint32_t>(count);
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
        const auto command=TraceCommandLine::parse(argc,argv);
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
        tracer.close();
        std::cout<<"Saved "<<command.output_path<<" (outgoing vertices, NOT a phase-function LUT).\n";
        if(errors){std::cerr<<"Unresolved vertices remain; inspect status/diagnostics before patch construction.\n";return EXIT_FAILURE;}
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
