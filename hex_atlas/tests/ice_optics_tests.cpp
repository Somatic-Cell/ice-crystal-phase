#include <ice_crystal/ice_optics.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <stdexcept>
using iceCrystal::optics::evaluate_ice_wb2008;
void check(bool value){if(!value)throw std::runtime_error("ice optical test failed");}
int main(int argc, char** argv)
{
    try {
        check(evaluate_ice_wb2008(450).n==1.3157);
        check(evaluate_ice_wb2008(550).n==1.3110);
        check(evaluate_ice_wb2008(650).n==1.3080);
        const auto s=evaluate_ice_wb2008(std::sqrt(550.0*560.0));
        check(std::abs(s.n-(1.3110+1.3106)/2)<4e-16);
        check(std::abs(s.kappa/std::sqrt(2.289e-9*2.839e-9)-1)<5e-15);
        check(evaluate_ice_wb2008(200).lower_nm==199);
        check(evaluate_ice_wb2008(1100).n==1.2998);
        for(double x:{199.,1101.,0.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
            bool caught=false;try{(void)evaluate_ice_wb2008(x);}catch(const std::invalid_argument&){caught=true;}check(caught);
        }
        double last=evaluate_ice_wb2008(200).n;
        for(int w=201;w<=1100;++w){const auto p=evaluate_ice_wb2008(w);check(p.n<last);last=p.n;}
        if(argc==2) {
            std::ifstream input(argv[1]);check(bool(input));
            std::string line;unsigned count=0;
            while(std::getline(input,line)) {
                if(line.empty()||line.front()=='#')continue;
                std::istringstream row(line);double um=0,n=0,k=0;check(bool(row>>um>>n>>k));
                const double nm=1000*um;
                if(nm<200||nm>1100)continue;
                const auto v=evaluate_ice_wb2008(nm);
                check(std::abs(v.n-n)<1e-15);check(std::abs(v.kappa/k-1)<5e-13);++count;
            }
            check(count==78);
        }
        std::cout<<"WB2008 exact knots / log-wavelength interpolation / domain: PASS\n";
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
