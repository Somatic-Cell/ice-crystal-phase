#include <rainbow/water_refractive_index.hpp>

#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>

// Coefficients/formulations: International Association for the Properties of
// Water and Steam (IAPWS), R9-97 (1997), Table 1 and Eq. (1); R6-95(2018),
// Tables 2-3 and Eq. (6). Both releases permit reproduction with attribution.
// https://iapws.org/documents/release/Rindex.download
// https://iapws.org/documents/release/IAPWS-95.download
// This code evaluates ALL 56 residual terms, not a polynomial density surrogate.
namespace rainbow
{
namespace
{
struct Sum
{
    double value = 0.0, correction = 0.0;
    void add(double x) noexcept
    {
        const double t = value + x;
        correction += std::abs(value) >= std::abs(x)
            ? (value-t)+x : (x-t)+value;
        value=t;
    }
    [[nodiscard]] double total() const noexcept { return value+correction; }
};
struct Term { double n, t; unsigned d, c; }; // c=0 means no exponential
constexpr std::array<Term, 51> terms{{
    { .12533547935523e-1, -.5, 1,0},
    { .78957634722828e1, .875,1,0},
    {-.87803203303561e1, 1,1,0},
    { .31802509345418, .5,2,0},
    {-.26145533859358, .75,2,0},
    {-.78199751687981e-2, .375,3,0},
    { .88089493102134e-2, 1,4,0},
    {-.66856572307965, 4,1,1},
    { .20433810950965, 6,1,1},
    {-.66212605039687e-4, 12,1,1},
    {-.19232721156002, 1,2,1},
    {-.25709043003438, 5,2,1},
    { .16074868486251, 4,3,1},
    {-.40092828925807e-1, 2,4,1},
    { .39343422603254e-6, 13,4,1},
    {-.75941377088144e-5, 9,5,1},
    { .56250979351888e-3, 3,7,1},
    {-.15608652257135e-4, 4,9,1},
    { .11537996422951e-8, 11,10,1},
    { .36582165144204e-6, 4,11,1},
    {-.13251180074668e-11, 13,13,1},
    {-.62639586912454e-9, 1,15,1},
    {-.10793600908932, 7,1,2},
    { .17611491008752e-1, 1,2,2},
    { .22132295167546, 9,2,2},
    {-.40247669763528, 10,2,2},
    { .58083399985759, 10,3,2},
    { .49969146990806e-2, 3,4,2},
    {-.31358700712549e-1, 7,4,2},
    {-.74315929710341, 10,4,2},
    { .47807329915480, 10,5,2},
    { .20527940895948e-1, 6,6,2},
    {-.13636435110343, 10,6,2},
    { .14180634400617e-1, 10,7,2},
    { .83326504880713e-2, 1,9,2},
    {-.29052336009585e-1, 2,9,2},
    { .38615085574206e-1, 3,9,2},
    {-.20393486513704e-1, 4,9,2},
    {-.16554050063734e-2, 8,9,2},
    { .19955571979541e-2, 6,10,2},
    { .15870308324157e-3, 9,10,2},
    {-.16388568342530e-4, 8,12,2},
    { .43613615723811e-1, 16,3,3},
    { .34994005463765e-1, 22,4,3},
    {-.76788197844621e-1, 23,4,3},
    { .22446277332006e-1, 23,5,3},
    {-.62689710414685e-4, 10,14,4},
    {-.55711118565645e-9, 50,3,6},
    {-.19905718354408, 44,6,6},
    { .31777497330738, 46,6,6},
    {-.11841182425981, 50,6,6}
}};

// Computes delta*d(phi^r)/d(delta). Using this directly avoids dividing each
// residual term by delta and multiplying the complete sum by delta again.
[[nodiscard]] double residual_pressure_factor(double delta, double tau)
{
    Sum sum{};
    for(const auto& a : terms)
    {
        const double dc = a.c ? std::pow(delta, double(a.c)) : 0.0;
        const double term = a.n * std::pow(delta, double(a.d)) * std::pow(tau,a.t)
                          * (a.c ? std::exp(-dc) : 1.0);
        sum.add(term * (double(a.d)-double(a.c)*dc));
    }
    struct Gaussian { double n,t,alpha,beta,gamma; };
    constexpr std::array<Gaussian,3> gaussians{{
        {-.31306260323435e2,0,20,150,1.21},
        { .31546140237781e2,1,20,150,1.21},
        {-.25213154341695e4,4,20,250,1.25}
    }};
    for(const auto& a : gaussians)
    {
        const double x=delta-1.0, y=tau-a.gamma;
        const double term=a.n*delta*delta*delta*std::pow(tau,a.t)
                         *std::exp(-a.alpha*x*x-a.beta*y*y);
        sum.add(term*(3.0-2.0*a.alpha*delta*x));
    }
    struct Critical { double n,b,C,D; };
    constexpr std::array<Critical,2> critical{{
        {-.14874640856724,.85,28,700},
        { .31806110878444,.95,32,800}
    }};
    // This private function is called only on the high-density liquid bracket,
    // far from delta=1. No critical-point special-case approximation is needed.
    const double x=delta-1.0, x2=x*x;
    constexpr double A=.32, B=.2, a=3.5, beta=.3;
    const double theta=1.0-tau+A*std::pow(x2,1.0/(2.0*beta));
    const double theta_d=A/beta*x*std::pow(x2,1.0/(2.0*beta)-1.0);
    const double Delta=theta*theta+B*std::pow(x2,a);
    const double Delta_d=2.0*theta*theta_d+2.0*B*a*x*std::pow(x2,a-1.0);
    for(const auto& v : critical)
    {
        const double psi=std::exp(-v.C*x2-v.D*(tau-1.0)*(tau-1.0));
        const double term=v.n*delta*std::pow(Delta,v.b)*psi;
        sum.add(term*(1.0+delta*(v.b*Delta_d/Delta-2.0*v.C*x)));
    }
    return sum.total();
}
[[nodiscard]] double pressure(double T, double rho)
{
    // The release's R, NOT an updated CODATA gas constant, is mandatory here.
    constexpr double R=461.51805, Tc=647.096, rhoc=322.0;
    return rho*R*T*(1.0+residual_pressure_factor(rho/rhoc,Tc/T));
}
}

double water_refractive_index_r9_97(double lambda_nm, double T, double rho)
{
    if(!(std::isfinite(lambda_nm) && lambda_nm>=200.0 && lambda_nm<=1100.0
        && std::isfinite(T) && T>=261.15 && T<=773.15
        && std::isfinite(rho) && rho>=0.0 && rho<=1060.0))
        throw std::invalid_argument("IAPWS R9-97: wavelength/state outside supported range.");
    const double l=lambda_nm/589.0, l2=l*l, t=T/273.15, r=rho/1000.0;
    constexpr double uv=.2292020, ir=5.432937;
    const double a=r*(.244257733+.00974634476*r-.00373234996*t
        +.000268678472*l2*t+.00158920570/l2
        +.00245934259/(l2-uv*uv)+.900704920/(l2-ir*ir)
        -.0166626219*r*r);
    const double n2=(1.0+2.0*a)/(1.0-a);
    if(!(n2>0.0 && std::isfinite(n2)))
        throw std::runtime_error("IAPWS R9-97 returned invalid index.");
    return std::sqrt(n2);
}

WaterOpticalState evaluate_water_optics(double lambda_nm, double celsius, double p)
{
    if(!(std::isfinite(celsius) && celsius>=0.0 && celsius<=60.0
        && std::isfinite(p) && p>=50000.0 && p<=100000000.0))
        throw std::invalid_argument("Water material supports liquid branch at 0..60 C and 50 kPa..100 MPa.");
    // Validate wavelength before solving the density equation.
    (void)water_refractive_index_r9_97(lambda_nm,celsius+273.15,1000.0);
    const double T=celsius+273.15;
    double lo=970.0, hi=1060.0;
    double flo=pressure(T,lo)-p, fhi=pressure(T,hi)-p;
    if(!(std::isfinite(flo) && std::isfinite(fhi) && flo<0.0 && fhi>0.0))
        throw std::runtime_error("IAPWS-95: liquid density bracket unavailable.");
    // Bisection in a predefined high-density branch. A pressure-only test is not
    // used as a convergence criterion; both bracket width and residual are kept.
    for(unsigned iteration=0; iteration<100; ++iteration)
    {
        const double mid=std::midpoint(lo,hi), fm=pressure(T,mid)-p;
        if(!std::isfinite(fm))throw std::runtime_error("IAPWS-95: nonfinite pressure.");
        if(fm<0.0){lo=mid;flo=fm;}else{hi=mid;fhi=fm;}
        if(hi-lo<=2e-10)
        {
            const bool use_lo=std::abs(flo)<std::abs(fhi);
            const double rho=use_lo?lo:hi, residual=use_lo?flo:fhi;
            if(std::abs(residual)>0.02)
                throw std::runtime_error("IAPWS-95: density pressure residual is too large.");
            return {T,p,rho,residual,hi-lo,water_refractive_index_r9_97(lambda_nm,T,rho)};
        }
    }
    throw std::runtime_error("IAPWS-95: liquid density solve did not converge.");
}
} // namespace rainbow
