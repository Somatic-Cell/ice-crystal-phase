#include <rainbow/phase_cdf.hpp>
#include <rainbow/phase_cdf_math.hpp>
#include <rainbow/phase_storage_math.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
using namespace rainbow;
namespace sm=phase_storage_math;
void require(bool b,const char* s){if(!b)throw std::runtime_error(s);}
void run(CudaContext& context,PhaseCdf& cdf,double sigma_degrees,bool coarse)
{
    constexpr unsigned nt=32,np=64;
    const unsigned ot=coarse?8:nt,op=coarse?16:np;
    const auto e=make_phase_u_edges(nt);
    std::vector<double> d(nt*np),fine(nt*np),filtered(nt*np),coarse_d(ot*op),ce(ot+1);
    for(unsigned i=0;i<nt;++i)for(unsigned j=0;j<np;++j)
        d[i*np+j]=.3+3*(1-e[i]-e[i+1])*(1-e[i]-e[i+1])+.05*i+.1*(j%7);
    const double maximum=*std::max_element(d.begin(),d.end());
    for(std::size_t k=0;k<d.size();++k)fine[k]=d[k]/maximum;
    sm::Grid grid{fine.data(),e.data(),nt,np};
    for(unsigned i=0;i<nt;++i)for(unsigned j=0;j<np;++j)
        filtered[i*np+j]=sigma_degrees>0?sm::gaussian_at(grid,i,j,sigma_degrees*sm::pi/180.,4):fine[i*np+j];
    grid.values=filtered.data();
    for(unsigned i=0;i<=ot;++i)ce[i]=e[i*(nt/ot)];
    for(unsigned i=0;i<ot;++i)for(unsigned j=0;j<op;++j)
        coarse_d[i*op+j]=sm::coarsen_at(grid,i,j,ot,op);
    DeviceBuffer<double> gpu(context);gpu.allocate(d.size());
    RAINBOW_CUDA_CHECK(cuMemcpyHtoD(gpu.address(),d.data(),gpu.byte_size()));
    PhaseDensityView v{};v.scalar=gpu.data();v.count=d.size();
    PhaseStorageSettings settings{};settings.theta_count=ot;settings.phi_count=op;
    settings.gaussian_sigma_degrees=sigma_degrees;
    // This test isolates GPU/CPU agreement; it does not certify optical resolution.
    PhaseCdfPolicy policy{};policy.allow_underresolved=true;
    cdf.build(v,nt,np,policy,settings);
    require(cdf.valid()&&cdf.theta_count()==ot&&cdf.phi_count()==op,"Wrong saved grid");
    std::vector<double> a(cdf.phi_cdf().element_count()),b(cdf.theta_cdf().element_count());
    cdf.phi_cdf().download(a);cdf.theta_cdf().download(b);
    long double total=0,moment=0;double l1=0;
    for(unsigned i=0;i<ot;++i)for(unsigned j=0;j<op;++j)
        total+=static_cast<long double>(coarse_d[i*op+j])*(ce[i+1]-ce[i]);
    for(unsigned i=0;i<ot;++i)for(unsigned j=0;j<op;++j)
    {
        const auto k=std::size_t(j)*(ot+1u)+i;
        const double recovered=(a[j+1]-a[j])*(b[k+1]-b[k]);
        const auto expected=static_cast<long double>(coarse_d[i*op+j])*(ce[i+1]-ce[i])/total;
        l1+=std::abs(recovered-static_cast<double>(expected));
        moment+=expected*sm::mean_cosine(ce[i],ce[i+1]);
    }
    require(l1<2e-10,"GPU storage CDF differs from CPU pipeline");
    require(std::abs(cdf.storage_statistics().g_stored-static_cast<double>(moment))<2e-10,"GPU HG label mismatch");
    require(std::abs(cdf.storage_statistics().aggregation_integral_relative_change)<1e-10,"Aggregation changed mass");
    require(std::isfinite(cdf.storage_statistics().coarsening_tv),"Missing coarsening error");
    if(!coarse)require(cdf.storage_statistics().coarsening_tv==0,"False coarsening error");
    // Invalid data MUST be rejected before smoothing, not concealed by it.
    d[37]=std::numeric_limits<double>::quiet_NaN();
    RAINBOW_CUDA_CHECK(cuMemcpyHtoD(gpu.address(),d.data(),gpu.byte_size()));
    bool rejected=false;try{cdf.build(v,nt,np,policy,settings);}catch(const std::runtime_error&){rejected=true;}
    require(rejected&&!cdf.valid(),"Prefilter concealed NaN");
    std::cout<<"Storage GPU sigma="<<sigma_degrees<<" coarse="<<coarse<<" L1="<<l1<<" passed\n";
}
}
#if defined(_WIN32)
int wmain(int argc,wchar_t** argv)
#else
int main(int argc,char** argv)
#endif
{
    try
    {
        if(argc!=2)throw std::invalid_argument("Expected phase_cdf.fatbin");
        CudaContext context;PhaseCdf cdf(context);cdf.load_module(std::filesystem::path(argv[1]));
        run(context,cdf,0,false);run(context,cdf,0,true);
        run(context,cdf,12,false);run(context,cdf,12,true);
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
