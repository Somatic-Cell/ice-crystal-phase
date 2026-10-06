#include <rainbow/patch_failure_report.hpp>
#include <rainbow/raindrop_tracer.hpp>
#include <rainbow/patch_accel.hpp>
#include <rainbow/patch_query.hpp>
#include <rainbow/patch_optics.hpp>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <span>
#include <stdexcept>
#include <vector>

namespace
{
void require(bool condition,const char* message)
{if(!condition)throw std::runtime_error(message);}
template<class T> std::vector<T> download(const rainbow::DeviceBuffer<T>& buffer)
{std::vector<T> out(buffer.element_count());buffer.download(std::span<T>{out});return out;}
template<class T> bool same(const std::vector<T>& a,const std::vector<T>& b)
{return a.size()==b.size()&&(a.empty()||std::memcmp(a.data(),b.data(),a.size()*sizeof(T))==0);}
}
int main(int argc,char** argv)
{
    try
    {
        if(argc!=5)throw std::invalid_argument("Expected raindrop.optixir patch_build.fatbin patch_query.optixir patch_optics.fatbin");
        rainbow::CudaContext cuda{0};rainbow::RaindropTracer tracer(cuda);
        rainbow::RaindropSettings settings{};
        settings.radius_mm=1.0f;settings.grid_width=settings.grid_height=129;
        const float inclination=20.0f*0.01745329251994329577f;
        settings.incident_direction={std::cos(inclination),-std::sin(inclination),0};
        tracer.create_pipeline(std::filesystem::path(argv[1]));tracer.trace(settings);
        rainbow::PatchAccel patches(cuda,tracer.optix_context());
        patches.load_module(std::filesystem::path(argv[2]));patches.build(tracer.vertices(),tracer.config());
        rainbow::PatchQuery query(cuda,tracer.optix_context());
        query.create_pipeline(std::filesystem::path(argv[3]));query.query_grid(patches,tracer.config(),90,180);
        rainbow::PatchOptics optics(cuda);optics.load_module(std::filesystem::path(argv[4]));
        rainbow::WaveOpticsSettings wave{};optics.evaluate_wave(patches,query,tracer.config(),wave);
        const auto before_vertices=download(tracer.vertices());
        const auto before_patches=download(patches.patches());
        const auto before_status=download(patches.cell_statuses());
        const auto before_hits=download(query.hits());const auto before_optics=download(optics.results());
        const auto report=rainbow::PatchFailureReport::capture(cuda,patches,query,optics,tracer.config(),wave);
        require(same(before_vertices,download(tracer.vertices())),"Capture altered source vertices");
        require(same(before_patches,download(patches.patches())),"Capture altered patch records");
        require(same(before_status,download(patches.cell_statuses())),"Capture altered cell classification");
        require(same(before_hits,download(query.hits())),"Capture altered query hits");
        require(same(before_optics,download(optics.results())),"Capture altered optical outputs");
        require(report.snapshot().origin==rainbow::PatchWitnessOrigin::GpuCapture,"Wrong provenance");
        require(report.snapshot().focal_unavailable_directions>=optics.wave_statistics().focal_errors,"Focal summary mismatch");
        // No exact count of pending patches is assumed: future correct classifier
        // improvements may recover them. Capture must still work when none fail.
        const auto path=std::filesystem::temp_directory_path()/
            ("rainbow_patch_capture_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".json");
        report.write_json(path);require(std::filesystem::file_size(path)>0,"Empty report");
        std::filesystem::remove(path);
        optics.close();query.close();patches.close();tracer.close();
        std::cout<<"Sparse GPU witness capture preserved all checked device buffers.\n";
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
