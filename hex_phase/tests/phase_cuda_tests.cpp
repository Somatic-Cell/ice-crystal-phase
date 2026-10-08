#include <ice_crystal/phase_cuda.hpp>
#include <ice_crystal/phase_generate.hpp>
#include <rainbow/phase_cdf.hpp>
#include <rainbow/cuda_error.hpp>
#include <rainbow/cuda_host_upload.hpp>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>

using namespace iceCrystal;
namespace
{
void check(bool b,const char* m){if(!b)throw std::runtime_error(m);}
void near(double x,double y,double e,const char* m)
{check(finite_value(x)&&finite_value(y)&&std::fabs(x-y)<=e,m);}
void compare_cdf(rainbow::CudaContext& context,const std::filesystem::path& modules,
    const PhaseGrid& grid,const std::vector<double>& mass,const rainbow::DeviceBuffer<double>& density,
    unsigned out_nt,unsigned out_np)
{
    std::cout << "  [cdf-case] " << grid.nt << 'x' << grid.np
              << " -> " << out_nt << 'x' << out_np << std::endl;
    const auto cpu=build_host_cdf(grid,mass,out_nt,out_np);
    rainbow::PhaseCdf gpu(context);gpu.load_module(modules/"phase_cdf.fatbin");
    rainbow::PhaseDensityView input{};input.count=density.element_count();
    input.scalar=reinterpret_cast<const double*>(density.address());input.stage=rainbow::PhaseDensityStage::Scalar;
    rainbow::PhaseStorageSettings settings{};settings.theta_count=out_nt;settings.phi_count=out_np;
    settings.gaussian_sigma_degrees=0;
    gpu.build(input,grid.nt,grid.np,{},settings);
    check(gpu.u_edges()==cpu.grid.edges,"CPU/CUDA saved u edges differ");
    std::vector<double> phi(out_np+1),theta(std::size_t(out_np)*(out_nt+1));
    gpu.phi_cdf().download(phi);gpu.theta_cdf().download(theta);
    Sum l1;
    for(unsigned j=0;j<out_np;++j)for(unsigned i=0;i<out_nt;++i)
    {
        const auto k=std::size_t(j)*(out_nt+1ull)+i;
        const double a=(phi[j+1]-phi[j])*(theta[k+1]-theta[k]);
        const double b=(cpu.phi[j+1]-cpu.phi[j])*(cpu.theta[k+1]-cpu.theta[k]);
        l1.add(std::fabs(a-b));
    }
    check(l1.value<2e-10,"CPU/CUDA cell probability L1 mismatch");
    near(gpu.storage_statistics().g_stored,cpu.report.g_stored,2e-11,"CPU/CUDA HG mismatch");
    near(gpu.storage_statistics().coarsening_tv,cpu.report.coarsening_tv,2e-10,"CPU/CUDA coarsening TV mismatch");
    std::cout<<"  CDF_L1="<<l1.value<<" g="<<gpu.storage_statistics().g_stored<<'\n';
}
// Bypass the hexagonal tracer and histogram so CDF storage errors are isolated.
// Repeat allocation/upload/coarsening on a non-blocking stream for both grids.
void direct_cdf(rainbow::CudaContext& context,const std::filesystem::path& modules)
{
    for(unsigned repeat=0;repeat<3;++repeat) for(unsigned nt:{72u,74u})
    {
        const unsigned np=nt==74 ? 146u : 144u;
        const auto grid=make_phase_grid(nt,np,{1,-0.4,0.3});
        for(bool poles:{false,true})
        {
            std::cout << "[direct-cdf] repeat=" << repeat << " poles=" << poles << std::endl;
            std::vector<double> mass(std::size_t(nt)*np);
            for(unsigned i=0;i<nt;++i) for(unsigned j=0;j<np;++j)
                if(!poles || i==0 || i+1==nt)
                    mass[std::size_t(i)*np+j]=(4*pi/np)*(grid.edges[i+1]-grid.edges[i]);
            const auto values=mass_to_density(grid,mass);
            rainbow::DeviceBuffer<double> device(context);device.allocate(values.size());
            rainbow::upload_host_to_device_sync(context,device.address(),values.data(),device.byte_size());
            compare_cdf(context,modules,grid,mass,device,nt,np);
            compare_cdf(context,modules,grid,mass,device,nt/2,np/2);
        }
    }
}
void synthetic(rainbow::CudaContext& context,const std::filesystem::path& modules)
{
    const auto grid=make_phase_grid(74,146,{1,-0.4,0.3});
    for(bool only_poles:{false,true})
    {
        std::cout << "[histogram] synthetic only_poles=" << only_poles << std::endl;
        PhaseMassAccumulator cpu(grid);PhaseMassAccumulatorCuda gpu(context,grid);
        gpu.load_module(modules/"phase_accumulate.fatbin");
        std::vector<OutgoingSample> samples(8192);
        for(std::size_t i=0;i<samples.size();++i)
        {
            auto& s=samples[i];
            if(only_poles||i%17==0)s.direction=(i%3)?grid.ki:-grid.ki;
            else
            {
                const double u=sample_uniform(7,i,0),phi=2*pi*sample_uniform(7,i,1)-pi;
                const double st=2*std::sqrt(u*(1-u));
                s.direction=grid.ki*(1-2*u)+grid.frame.e0*(st*std::cos(phi))+grid.frame.e1*(st*std::sin(phi));
            }
            s.power_fraction=std::pow(10.,-12*sample_uniform(7,i,2));
        }
        rainbow::DeviceBuffer<OutgoingSample> device(context);device.allocate(samples.size());context.make_current();
        rainbow::upload_host_to_device_sync(context, device.address(), samples.data(), device.byte_size());
        cpu.add(samples,{},0.7);gpu.add(device,{},0.7);gpu.finish();
        const auto mass=cpu.masses();const auto d=gpu.download_density();const auto expected=mass_to_density(grid,mass);
        Sum error;for(std::size_t k=0;k<d.size();++k)
        {const auto i=k/grid.np;error.add(std::fabs(d[k]-expected[k])*(4*pi/grid.np)*(grid.edges[i+1]-grid.edges[i]));}
        check(error.value/cpu.statistics().point_mass<2e-11,"GPU histogram differs from compensated CPU");
        check(gpu.statistics().samples==samples.size(),"lost GPU samples");
        compare_cdf(context,modules,grid,mass,gpu.density(),37,73);
    }
    // A bad sample cannot become a valid empty/sanitized histogram.
    PhaseMassAccumulatorCuda broken(context,grid);broken.load_module(modules/"phase_accumulate.fatbin");
    OutgoingSample invalid{};invalid.direction={0,0,0};invalid.power_fraction=1;
    rainbow::DeviceBuffer<OutgoingSample> d(context);d.allocate(1);context.make_current();
    rainbow::upload_host_to_device_sync(context, d.address(), &invalid, sizeof(invalid));broken.add(d,{},1);
    bool rejected=false;try{broken.finish();}catch(const std::exception&){rejected=true;}
    check(rejected,"GPU silently accepted invalid sample");
}
void real(rainbow::CudaContext& context,const std::filesystem::path& modules)
{
    std::cout << "[trace] weighted orientations" << std::endl;
    const auto plan=parse_orientation_plan("qw,qx,qy,qz,weight,samples,seed\n1,0,0,0,0.25,256,7\n"
        "0.7071067811865476,0,0,0.7071067811865476,0.75,512,13\n");
    PhaseGenerationSettings s;s.prism={.1,.2};s.grid=make_phase_grid(74,146,{1,-.4,.3});
    s.interior_index=1.31;s.wavelength_nm=550;s.out_nt=37;s.out_np=73;
    PhaseMassAccumulator cpu(s.grid);PhaseMassAccumulatorCuda gpu(context,s.grid);
    gpu.load_module(modules/"phase_accumulate.fatbin");HexPrismTracer tracer(context);tracer.load_module(modules/"hex_trace.fatbin");
    CompositionAudit audit;
    for(const auto& n:plan.nodes)
    {
        const auto c=node_trace_settings(s,n);const double area=node_area(s,n),scale=n.weight*area/double(n.samples);
        audit.expected_area.add(n.weight*area);
        for(std::uint64_t first=0;first<n.samples;)
        {
            const auto count=std::min<std::uint64_t>(73,n.samples-first);
            const auto cb=trace_batch_cpu(c,first,count);require_accepted(summarize(cb.audits));cpu.add(cb.outgoing,n.rotation,scale);
            tracer.trace(c,first,count);record_batch_audit(audit,summarize_phase_rays(tracer.download_audits()),scale);
            gpu.add(tracer.outgoing(),n.rotation,scale);first+=count;
        }
        ++audit.orientations_done;
    }
    gpu.finish();audit.histogram=gpu.statistics();audit.grid_moments=gpu.moments();validate_composition(audit,s,plan);
    compare_cdf(context,modules,s.grid,cpu.masses(),gpu.density(),37,73);
    std::cout<<"  real_histogram_mass_error="<<audit.histogram_relative_error<<'\n';
}
}
int main(int argc,char** argv)
{
    try
    {
        const auto modules=argc==2?std::filesystem::path(argv[1]):std::filesystem::absolute(argv[0]).parent_path()/"modules";
        std::cout<<std::setprecision(17);rainbow::CudaContext context(0);
        direct_cdf(context,modules);synthetic(context,modules);real(context,modules);
        std::cout<<"PASS CUDA histogram + existing PhaseCdf: synthetic HDR, poles, invalid sample, real weighted orientations\n";return 0;
    }
    catch(const std::exception& e){std::cerr<<"CUDA test FAILED (not skipped): "<<e.what()<<'\n';return 1;}
}
