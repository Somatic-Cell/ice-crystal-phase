#include <rainbow/projected_area.hpp>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{
unsigned checks=0;
void require(bool b,const char* message){++checks;if(!b)throw std::runtime_error(message);}
void near(double a,double b,double rel=2e-10){require(std::abs(a-b)<=rel*std::abs(b)+1e-13,"Area mismatch");}
const std::array<std::array<float,8>,6> rows{{
    {{0,0,0,0,0,0,0,0}},
    {{-.0131f,-.0120f,-.0376f,-.0096f,-.0004f,.0015f,.0005f,0}},
    {{-.0282f,-.0230f,-.0779f,-.0175f,.0021f,.0046f,.0011f,-.0006f}},
    {{-.0458f,-.0335f,-.1211f,-.0227f,.0083f,.0089f,.0012f,-.0021f}},
    {{-.0644f,-.0416f,-.1629f,-.0246f,.0176f,.0131f,.0002f,-.0044f}},
    {{-.0840f,-.0480f,-.2034f,-.0237f,.0297f,.0166f,-.0021f,-.0072f}}
}};
}
int main(int argc,char** argv)
{
    try
    {
        // Small independently callable probe for Python projection references.
        if(argc==3)
        {
            const int row=std::stoi(argv[1]);const double angle=std::stod(argv[2]);
            if(row<0||row>=6)throw std::invalid_argument("row");
            std::array<double,8> c{};for(unsigned i=0;i<8;++i)c[i]=rows[static_cast<unsigned>(row)][i];
            const auto r=rainbow::compute_projected_area(c,1,{std::cos(angle),-std::sin(angle),0});
            std::cout.precision(17);std::cout<<r.area_mm2<<' '<<r.estimated_absolute_error_mm2<<' '<<r.convexity_intervals<<' '<<r.quadrature_evaluations<<'\n';return 0;
        }
        constexpr double pi=3.141592653589793238462643383279502884;
        for(double a:{.1,.4,1.0,3.0})
        for(auto k:{std::array<double,3>{1,0,0},{0,-1,0},{1,2,3}})
        {
            const auto r=rainbow::compute_projected_area({},a,k);near(r.area_mm2,pi*a*a,1e-14);
            require(r.analytic_sphere,"Sphere not analytic");
        }
        // All tabulated profiles AND interpolated profiles; coefficient values
        // are binary32 as in RaindropShape, not a different double-precision model.
        for(unsigned j=0;j+1<rows.size();++j)for(unsigned step=0;step<=10;++step)
        {
            const float t=float(step)/10;
            std::array<double,8> c{};
            for(unsigned n=0;n<8;++n)c[n]=std::fma(t,rows[j+1][n]-rows[j][n],rows[j][n]);
            for(double angle:{-.9,0.,.3490658503988659,1.5707963267948966})
            {
                const std::array<double,3> k{std::cos(angle),-std::sin(angle),0};
                const auto r=rainbow::compute_projected_area(c,1,k);
                require(r.area_mm2>0 && r.estimated_relative_error<1.01e-10,"Invalid area or error");
                near(r.area_mm2,rainbow::compute_projected_area(c,1,{-k[0],-k[1],-k[2]}).area_mm2,1e-14);
                near(r.area_mm2,rainbow::compute_projected_area(c,1,{0,k[1],k[0]}).area_mm2,1e-14);
                near(4*r.area_mm2,rainbow::compute_projected_area(c,2,k).area_mm2,1e-14);
            }
        }
        for(int kind=0;kind<4;++kind)
        {
            bool rejected=false;
            try
            {
                std::array<double,8> c{};double a=1;std::array<double,3> k{1,0,0};
                if(kind==0)c[2]=.8; // nonconvex radial profile
                if(kind==1)c[0]=-2;
                if(kind==2)k={0,0,0};
                if(kind==3)c[1]=std::numeric_limits<double>::quiet_NaN();
                (void)rainbow::compute_projected_area(c,a,k);
            }catch(const std::exception&){rejected=true;}
            require(rejected,"Invalid geometry accepted");
        }
        std::cout<<"projected_area checks="<<checks<<" PASS\n";
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
