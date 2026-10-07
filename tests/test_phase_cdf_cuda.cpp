#include <rainbow/phase_cdf.hpp>
#include <rainbow/phase_cdf_math.hpp>
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
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
void run_case(CudaContext& context,PhaseCdf& cdf,unsigned nt,unsigned np)
{
    const auto edges=make_phase_u_edges(nt);
    std::vector<double> density(std::size_t(nt)*np);
    long double sum=0;
    for(unsigned i=0;i<nt;++i)for(unsigned j=0;j<np;++j)
    {
        double d=.2+(i+1.0)*(i+1.0)+.1*j;
        if(np>2 && (j==0 || j+1==np || (i+2*j)%5==0))d=0;
        density[std::size_t(i)*np+j]=d;
        sum+=static_cast<long double>(d)*static_cast<long double>(edges[i+1]-edges[i]);
    }
    DeviceBuffer<double> input(context);input.allocate(density.size());
    RAINBOW_CUDA_CHECK(cuMemcpyHtoD(input.address(),density.data(),input.byte_size()));
    PhaseDensityView view{};view.scalar=input.data();view.count=density.size();
    cdf.build(view,nt,np);
    require(cdf.valid(),"CDF was not published");
    std::vector<double> a(cdf.phi_cdf().element_count()),b(cdf.theta_cdf().element_count());
    cdf.phi_cdf().download(a);cdf.theta_cdf().download(b);
    require(a.front()==0&&a.back()==1,"Bad marginal endpoints");
    double l1=0;
    for(unsigned j=0;j<np;++j)for(unsigned i=0;i<nt;++i)
    {
        const auto k=std::size_t(j)*(nt+1)+i;
        const double recovered=(a[j+1]-a[j])*(b[k+1]-b[k]);
        const double expected=static_cast<double>(static_cast<long double>(density[std::size_t(i)*np+j])*
            static_cast<long double>(edges[i+1]-edges[i])/sum);
        require(recovered>=0,"Negative probability mass");
        l1+=std::abs(recovered-expected);
    }
    require(l1<2e-12,"GPU CDF differs from independent CPU mass calculation");
    long double moment=0;
    for(unsigned j=0;j<np;++j)for(unsigned i=0;i<nt;++i)
        moment+=static_cast<long double>((a[j+1]-a[j])*(b[j*(nt+1u)+i+1]-b[j*(nt+1u)+i]))*
            static_cast<long double>(1.0-(edges[i]+edges[i+1]));
    require(std::abs(cdf.storage_statistics().g_stored-static_cast<double>(moment))<5e-12,
        "Saved-CDF HG moment disagrees with CPU reconstruction");
    density[0]=-1;
    RAINBOW_CUDA_CHECK(cuMemcpyHtoD(input.address(),density.data(),input.byte_size()));
    bool rejected=false;
    try{cdf.build(view,nt,np);}catch(const std::runtime_error&){rejected=true;}
    require(rejected&&!cdf.valid(),"Invalid input became a valid CDF");
    std::cout<<"GPU CDF "<<nt<<'x'<<np<<" L1="<<l1<<" passed\n";
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
        run_case(context,cdf,1,1);run_case(context,cdf,17,13);
        run_case(context,cdf,257,33);run_case(context,cdf,513,513);
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
