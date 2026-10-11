#include <ice_crystal/atlas_cuda.hpp>
#include <algorithm>
#include <iostream>
#include <stdexcept>
int main(int argc,char** argv)
{
    (void)argc;
    try
    {
        rainbow::CudaContext context;
        const auto modules=std::filesystem::absolute(argv[0]).parent_path()/"modules";
        iceCrystal::atlas::CudaTracer gpu(context,modules/"atlas_trace.fatbin");
        iceCrystal::atlas::Settings s;s.prism={.1,.2};s.interior_index=1.31;s.wavelength_nm=550;
        s.incident_theta_intervals=3;s.incident_phi_intervals=2;s.theta_count=12;s.phi_count=24;s.area_subdivisions=2;
        for(unsigned i=0;i<=3;++i)for(unsigned j=0;j<=2;++j)
        {
            const auto n=iceCrystal::atlas::make_node(s,i,j);const auto cpu=iceCrystal::atlas::trace_cpu(s,n);
            s.batch_rays=17;const auto a=gpu.trace(s,n);const auto again=gpu.trace(s,n);
            if(a.mass_mm2!=again.mass_mm2)throw std::runtime_error("Repeated GPU atlas not bitwise equal");
            s.batch_rays=64;const auto b=gpu.trace(s,n);
            if(a.mass_mm2!=b.mass_mm2)throw std::runtime_error("GPU ordered sum depends on batch size");
            double error=0;for(std::size_t k=0;k<a.mass_mm2.size();++k)error+=std::fabs(a.mass_mm2[k]-cpu.mass_mm2[k]);
            if(error/cpu.projected_area_mm2>2e-9||a.audit.rays!=cpu.audit.rays||a.audit.outputs!=cpu.audit.outputs)
                throw std::runtime_error("CPU/CUDA atlas mass or trace count mismatch");
            std::cout<<"atlas_node="<<i<<','<<j<<" L1_relative="<<error/cpu.projected_area_mm2<<" repeat=bitwise batch=bitwise\n";
        }
        std::cout<<"CUDA atlas PASS (not skipped)\n";return 0;
    }
    catch(const std::exception& e){std::cerr<<"CUDA atlas FAILED (not skipped): "<<e.what()<<'\n';return 1;}
}
