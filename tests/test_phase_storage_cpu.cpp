#include <rainbow/phase_storage_math.hpp>
#include <rainbow/phase_cdf_math.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
using namespace rainbow;
namespace sm=phase_storage_math;
std::uint64_t checks=0;
void require(bool b,const char* s){++checks;if(!b)throw std::runtime_error(s);}
long double integral(const std::vector<double>& d,const std::vector<double>& e,unsigned nt,unsigned np)
{
    long double v=0;
    for(unsigned i=0;i<nt;++i)for(unsigned j=0;j<np;++j)
        v+=static_cast<long double>(d[std::size_t(i)*np+j])*(static_cast<long double>(e[i+1])-e[i])/np;
    return v;
}
long double reference_gaussian(sm::Grid g,unsigned row,unsigned col,double sigma,double support)
{
    // Independent all-pairs quadrature, no spherical-cap pruning / modulo window.
    const long double pi=acosl(-1.L), t=(row+.5L)*pi/g.nt, phi=2*pi*col/g.np;
    long double sum=0,w=0;
    for(unsigned r=0;r<g.nt;++r)for(unsigned c=0;c<g.np;++c)
    {
        const long double a=(r+.5L)*pi/g.nt,b=2*pi*c/g.np;
        const long double dx=sinl(t)*cosl(phi)-sinl(a)*cosl(b);
        const long double dy=sinl(t)*sinl(phi)-sinl(a)*sinl(b);
        const long double dz=cosl(t)-cosl(a);
        const long double angle=2*asinl(std::min(1.L,sqrtl(dx*dx+dy*dy+dz*dz)/2));
        if(angle>static_cast<long double>(sigma*support)+0x1p-44L)continue;
        const long double k=expl(-.5L*(angle/sigma)*(angle/sigma))*(static_cast<long double>(g.edges[r+1])-g.edges[r]);
        w+=k;sum+=k*g.values[std::uint64_t(r)*g.np+c];
    }
    return sum/w;
}
void test_gaussian()
{
    for(const auto dims:std::array<std::array<unsigned,2>,3>{{{9,18},{12,25},{17,1}}})
    {
        const unsigned nt=dims[0],np=dims[1];auto e=make_phase_u_edges(nt);
        std::vector<double> values(nt*np),out(nt*np);
        for(unsigned i=0;i<nt;++i)for(unsigned j=0;j<np;++j)
            values[std::size_t(i)*np+j]=.2+.03*i+.4*(j%3)+.02*(i*j%11);
        sm::Grid g{values.data(),e.data(),nt,np};const double sigma=16.*sm::pi/180.;
        for(unsigned i=0;i<nt;++i)for(unsigned j=0;j<np;++j)
        {
            const auto x=sm::gaussian_at(g,i,j,sigma,4);
            const auto y=reference_gaussian(g,i,j,sigma,4);
            require(std::isfinite(x)&&x>=0,"Gaussian produced invalid density");
            require(std::abs(static_cast<long double>(x)-y)<5e-12L,"Cap quadrature disagrees with independent reference");
            out[std::size_t(i)*np+j]=x;
        }
        // Cyclic longitude shifts must commute with the discrete operator.
        std::vector<double> rotated(values.size());
        for(unsigned i=0;i<nt;++i)for(unsigned j=0;j<np;++j)rotated[std::size_t(i)*np+j]=values[std::size_t(i)*np+(j+3)%np];
        sm::Grid rg{rotated.data(),e.data(),nt,np};
        for(unsigned i=0;i<nt;++i)for(unsigned j=0;j<np;++j)
            require(std::abs(sm::gaussian_at(rg,i,j,sigma,4)-out[std::size_t(i)*np+(j+3)%np])<1e-12,"Longitude seam error");
        std::fill(values.begin(),values.end(),2.7);
        for(unsigned i=0;i<nt;++i)for(unsigned j=0;j<np;++j)
            require(std::abs(sm::gaussian_at(g,i,j,sigma,4)-2.7)<3e-15,"Constant not preserved");
    }
}
void test_aggregation()
{
    for(const auto dims:std::array<std::array<unsigned,4>,4>{{{18,40,9,10},{32,64,8,16},{17,13,1,1},{8,1,4,1}}})
    {
        const auto nt=dims[0],np=dims[1],ot=dims[2],op=dims[3];
        auto e=make_phase_u_edges(nt);std::vector<double> d(nt*np),out(ot*op),oe(ot+1);
        for(unsigned i=0;i<=ot;++i)oe[i]=e[i*(nt/ot)];
        for(unsigned i=0;i<nt;++i)for(unsigned j=0;j<np;++j)d[i*np+j]=.1+(i+1.)*(i+1.)+.4*j;
        sm::Grid g{d.data(),e.data(),nt,np};
        for(unsigned i=0;i<ot;++i)for(unsigned j=0;j<op;++j)
        {
            const double x=sm::coarsen_at(g,i,j,ot,op);out[i*op+j]=x;
            long double m=0;
            for(unsigned r=i*(nt/ot);r<(i+1)*(nt/ot);++r)
                for(unsigned c=j*(np/op);c<(j+1)*(np/op);++c)
                    m+=static_cast<long double>(d[r*np+c])*(static_cast<long double>(e[r+1])-e[r]);
            m/=(static_cast<long double>(oe[i+1])-oe[i])*(np/op);
            require(std::abs(static_cast<long double>(x)-m)<2e-13L*std::max(1.L,m),"Cell mass mismatch");
        }
        const auto a=integral(d,e,nt,np),b=integral(out,oe,ot,op);
        require(std::abs(a-b)/a<3e-15,"Global mass not preserved by coarsening");
    }
}
void test_moments()
{
    for(unsigned nt:{1u,2u,17u,900u,1800u})
    {
        const auto e=make_phase_u_edges(nt);long double g=0;
        for(unsigned i=0;i<nt;++i)
        {
            const long double a=1-2*static_cast<long double>(e[i]);
            const long double b=1-2*static_cast<long double>(e[i+1]);
            const auto m=sm::mean_cosine(e[i],e[i+1]);
            require(std::abs(static_cast<long double>(m)-(a+b)/2)<4e-16L,"Wrong cell cosine average");
            g+=(static_cast<long double>(e[i+1])-e[i])*m;
        }
        require(std::abs(g)<5e-16L,"Uniform spherical distribution has nonzero g");
    }
    const double phi[3]={0,.25,1};
    const double c[2][4]={{0,.1,.5,1},{0,0,.25,1}};
    const double e[4]={0,.1,.7,1};long double m=0,g=0;
    for(unsigned j=0;j<2;++j)for(unsigned i=0;i<3;++i)
    {
        const double w=(phi[j+1]-phi[j])*(c[j][i+1]-c[j][i]);
        m+=w;g+=w*sm::mean_cosine(e[i],e[i+1]);
    }
    require(std::abs(m-1)<2e-16L,"CDF masses do not sum to one");
    require(std::abs(g-(-.40125L))<2e-16L,"Saved-CDF first moment mismatch");
}
}
int main()
{
    try{test_gaussian();test_aggregation();test_moments();std::cout<<checks<<" storage math checks passed.\n";return 0;}
    catch(const std::exception& e){std::cerr<<"Error after "<<checks<<" checks: "<<e.what()<<'\n';return 1;}
}
