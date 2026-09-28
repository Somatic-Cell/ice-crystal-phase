#include "raindrop_validation.hpp"
#include <rainbow/raindrop_tracer.hpp>
#include <cstdlib>

namespace
{
int validate(const std::filesystem::path& module)
{
    try
    {
        rainbow::CudaContext cuda_context{0};
        rainbow::RaindropTracer tracer(cuda_context);
        tracer.create_pipeline(module);
        for(const float radius:{0.4f,1.0f,3.0f})
        for(unsigned polarization=0;polarization<2;++polarization)
        {
            rainbow::RaindropSettings settings;
            settings.radius_mm=radius;settings.grid_width=settings.grid_height=17;
            settings.incident_direction={1,-0.3f,0.1f};
            settings.incident_field=polarization?rainbow::Field32{{0,0},{1,0}}:rainbow::Field32{{1,0},{0,0}};
            const auto config=settings.make_config();
            tracer.trace(settings); // 同じ module/pipeline/context を条件間で再利用する．
            const auto actual=tracer.download_vertices();
            if(radius<=0.4f)
            {
                rainbow::tests::AnalyticSphereIntersector reference;
                rainbow::tests::verify_trace(config,actual,rainbow::tests::cpu_trace(config,reference));
                rainbow::tests::verify_sphere_optical_lengths(config,actual);
            }
            else
            {
                rainbow::tests::CosineReferenceIntersector reference{config.shape};
                rainbow::tests::verify_trace(config,actual,rainbow::tests::cpu_trace(config,reference));
            }
        }
        tracer.close();
        std::cout<<"OptiX raindrop R/TT/TRT/TRRT generation: passed\n";
        return EXIT_SUCCESS;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return EXIT_FAILURE;}
}
}
#if defined(_WIN32)
int wmain(const int argc,wchar_t* argv[])
#else
int main(const int argc,char* argv[])
#endif
{
    if(argc!=2){std::cerr<<"Expected raindrop_trace.optixir path\n";return EXIT_FAILURE;}
    try{return validate(std::filesystem::path{argv[1]});}
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return EXIT_FAILURE;}
}
