#include "raindrop_validation.hpp"
#include <cstdlib>
int main()
{
    try
    {
        rainbow::tests::verify_shape_table();
        rainbow::tests::verify_vertical_orientation();
        rainbow::tests::verify_fresnel();
        rainbow::tests::verify_compensated_optical_path();
        for(float radius:{0.1f,0.4f,1.0f,2.0f,3.0f})
        {
            rainbow::RaindropSettings settings;
            settings.radius_mm=radius;settings.grid_width=settings.grid_height=17;
            settings.incident_direction={1.0f,-0.3f,0.1f};
            const auto config=settings.make_config();
            rainbow::RaindropIntersector solver{config.shape};
            const auto vertices=rainbow::tests::cpu_trace(config,solver);
            if(radius<=0.4f)
            {
                rainbow::tests::AnalyticSphereIntersector reference;
                rainbow::tests::verify_trace(config,vertices,rainbow::tests::cpu_trace(config,reference));
                rainbow::tests::verify_sphere_optical_lengths(config,vertices);
            }
            else
            {
                rainbow::tests::CosineReferenceIntersector reference{config.shape};
                rainbow::tests::verify_trace(config,vertices,rainbow::tests::cpu_trace(config,reference));
            }
        }
        std::cout<<"Raindrop CPU geometry/Fresnel/paths: passed\n";
        return EXIT_SUCCESS;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return EXIT_FAILURE;}
}
