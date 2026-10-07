#include <rainbow/projected_area.hpp>

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>

#ifdef __FAST_MATH__
#error projected_area.cpp requires IEEE arithmetic; do not use fast-math.
#endif

namespace rainbow
{
namespace
{
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double eps = std::numeric_limits<double>::epsilon();
static_assert(std::numeric_limits<double>::is_iec559 &&
              std::numeric_limits<double>::digits == 53);

struct Interval { double lo, hi; };
Interval point(double x) { return {x,x}; }
double down(double x) { return std::nextafter(x, -std::numeric_limits<double>::infinity()); }
double up(double x) { return std::nextafter(x, std::numeric_limits<double>::infinity()); }
Interval operator+(Interval x, Interval y) { return {down(x.lo+y.lo), up(x.hi+y.hi)}; }
Interval operator-(Interval x, Interval y) { return {down(x.lo-y.hi), up(x.hi-y.lo)}; }
Interval operator*(Interval x, Interval y)
{
    const std::array<double,4> v{x.lo*y.lo,x.lo*y.hi,x.hi*y.lo,x.hi*y.hi};
    return {down(*std::min_element(v.begin(),v.end())),
            up(*std::max_element(v.begin(),v.end()))};
}
Interval square(Interval x)
{
    const double hi = (std::max)(x.lo*x.lo,x.hi*x.hi);
    const double lo = x.lo<=0 && x.hi>=0 ? 0 : (std::min)(x.lo*x.lo,x.hi*x.hi);
    return {lo==0 ? 0 : down(lo), up(hi)};
}
bool finite(Interval x) { return std::isfinite(x.lo) && std::isfinite(x.hi); }

// Chebyshev recurrences and their first two derivatives in mu=cos(chi).
// Q=r^2+2(1-mu^2)r_mu^2-r[(1-mu^2)r_mumu-mu*r_mu] > 0
// and G=r+mu*r_mu > 0 are the two positive-curvature conditions.
// r>0 plus these inequalities gives a strictly convex radial surface of
// revolution. Pole limits use the polynomial expressions, without /sin(chi).
std::uint32_t check_convexity(const std::array<double,8>& c)
{
    struct Domain { double lo, hi; unsigned depth; };
    std::vector<Domain> stack{{-1,1,0}};
    std::uint32_t visited=0;
    while(!stack.empty())
    {
        if(++visited>262144u) throw std::runtime_error("Projected area: convexity verification work limit.");
        const auto d=stack.back(); stack.pop_back();
        const Interval mu{d.lo,d.hi};
        auto t0=point(1), t1=mu, p0=point(0), p1=point(1), s0=point(0), s1=point(0);
        auto r=point(1)+point(c[0])+point(c[1])*mu;
        auto rp=point(c[1]), rpp=point(0);
        for(unsigned n=2;n<8;++n)
        {
            const auto t=point(2)*mu*t1-t0;
            const auto p=point(2)*t1+point(2)*mu*p1-p0;
            const auto s=point(4)*p1+point(2)*mu*s1-s0;
            r=r+point(c[n])*t; rp=rp+point(c[n])*p; rpp=rpp+point(c[n])*s;
            t0=t1;t1=t;p0=p1;p1=p;s0=s1;s1=s;
        }
        auto sin2=point(1)-square(mu);
        sin2.lo=(std::max)(0.0,sin2.lo); // known range of 1-mu^2 on [-1,1]
        const auto g=r+mu*rp;
        const auto q=square(r)+point(2)*sin2*square(rp)-r*(sin2*rpp-mu*rp);
        if(!finite(r)||!finite(g)||!finite(q))
            throw std::runtime_error("Projected area: nonfinite convexity bound.");
        if(r.lo>0 && g.lo>0 && q.lo>0) continue;
        if(r.hi<=0 || g.hi<=0 || q.hi<=0)
            throw std::invalid_argument("Projected area: nonconvex or invalid radial shape; surface formula would double-count. No hull approximation is used.");
        if(d.depth==24)
            throw std::runtime_error("Projected area: could not verify strict convexity; no area was substituted.");
        const double mid=0.5*(d.lo+d.hi);
        stack.push_back({mid,d.hi,d.depth+1});
        stack.push_back({d.lo,mid,d.depth+1});
    }
    return visited;
}

struct Integrand
{
    const std::array<double,8>& c;
    double transverse, axial;
    double operator()(double chi) const
    {
        const double mu=std::cos(chi), sn=std::sin(chi);
        double t0=1,t1=mu,p0=0,p1=1;
        double r=1+c[0]+c[1]*mu, rp=c[1];
        for(unsigned n=2;n<8;++n)
        {
            const double t=2*mu*t1-t0, p=2*t1+2*mu*p1-p0;
            r+=c[n]*t;rp+=c[n]*p;t0=t1;t1=t;p0=p1;p1=p;
        }
        // n dA = r sin(chi) [(r sin- r_chi cos)e_rho +
        //                         (r cos+ r_chi sin)e_axis] dchi dphi.
        const double horizontal=sn*(r+mu*rp);
        const double vertical=mu*r-sn*sn*rp;
        const double a=std::abs(transverse*horizontal), b=std::abs(axial*vertical);
        // Exact azimuth integral J = integral_0^{2pi} |a*cos(phi)+b| dphi.
        double j;
        if(a<=b) j=2*pi*b;
        else j=4*(std::sqrt((a-b)*(a+b))+b*std::asin(b/a));
        const double value=0.5*r*sn*j;
        if(!(value>=0) || !std::isfinite(value))
            throw std::runtime_error("Projected area: invalid surface integrand.");
        return value;
    }
};

struct Panel
{
    double lo,hi,value,error;
    bool operator<(const Panel& other) const { return error<other.error; }
};
// Gauss(7)-Kronrod(15). The constants are quadrature nodes/weights; this is
// an independently written adaptive integrator. Error rescaling follows the
// usual embedded-rule variation/roundoff safeguards, not a rigorous bound.
Panel integrate_panel(const Integrand& f,double lo,double hi)
{
    constexpr std::array<double,8> x{
        .991455371120812639206854697526329,.949107912342758524526189684047851,
        .864864423359769072789712788640926,.741531185599394439863864773280788,
        .586087235467691130294144838258730,.405845151377397166906606412076961,
        .207784955007898467600689403773245,0};
    constexpr std::array<double,8> wk{
        .022935322010529224963732008058970,.063092092629978553290700663189204,
        .104790010322250183839876322541518,.140653259715525918745189590510238,
        .169004726639267902826583426598550,.190350578064785409913256402421014,
        .204432940075298892414161999234649,.209482141084727828012999174891714};
    constexpr std::array<double,4> wg{
        .129484966168869693270611432679082,.279705391489276667901467771423780,
        .381830050505118944950369775488975,.417959183673469387755102040816327};
    const double mid=0.5*(lo+hi), half=0.5*(hi-lo), fc=f(mid);
    double k=wk[7]*fc,g=wg[3]*fc;
    std::array<double,7> left{},right{};
    for(unsigned i=0;i<7;++i)
    {
        left[i]=f(mid-half*x[i]);right[i]=f(mid+half*x[i]);
        const double pair=left[i]+right[i]; k+=wk[i]*pair;
        if(i%2==1) g+=wg[i/2]*pair;
    }
    const double mean=0.5*k;
    double asc=wk[7]*std::abs(fc-mean);
    for(unsigned i=0;i<7;++i)asc+=wk[i]*(std::abs(left[i]-mean)+std::abs(right[i]-mean));
    asc*=half;
    const double value=k*half;
    double error=std::abs((k-g)*half);
    if(asc>0 && error>0) error=asc*(std::min)(1.0,std::pow(200*error/asc,1.5));
    error=(std::max)(error,50*eps*std::abs(value));
    return {lo,hi,value,error};
}
}

ProjectedAreaResult compute_projected_area(const std::array<double,8>& c,
    double radius,const std::array<double,3>& direction,const ProjectedAreaOptions& o)
{
    if(std::fegetround()!=FE_TONEAREST)
        throw std::runtime_error("Projected area requires round-to-nearest CPU arithmetic.");
    if(!(std::isfinite(radius)&&radius>0) ||
       !(std::isfinite(o.relative_tolerance)&&o.relative_tolerance>=1e-13&&o.relative_tolerance<1) ||
       !(std::isfinite(o.absolute_tolerance_dimensionless)&&o.absolute_tolerance_dimensionless>=0) ||
       o.maximum_quadrature_panels<32 || o.maximum_quadrature_panels>1000000)
        throw std::invalid_argument("Invalid projected-area options or radius.");
    for(double v:c)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite shape coefficient.");
    for(double v:direction)if(!std::isfinite(v))throw std::invalid_argument("Nonfinite incident direction.");
    const double norm=std::hypot(direction[0],direction[1],direction[2]);
    if(!(norm>0 && std::isfinite(norm)))throw std::invalid_argument("Invalid incident direction.");
    const double radius2=radius*radius;
    if(!(radius2>0 && std::isfinite(radius2)))throw std::invalid_argument("Radius squared is not representable.");
    ProjectedAreaResult result{};
    bool sphere=true;for(unsigned i=1;i<8;++i)sphere=sphere&&c[i]==0;
    if(sphere)
    {
        const double r=1+c[0];if(!(r>0&&std::isfinite(r)))throw std::invalid_argument("Invalid sphere radius.");
        result.area_mm2=pi*r*r*radius2;
        result.estimated_relative_error=8*eps;
        result.estimated_absolute_error_mm2=result.area_mm2*result.estimated_relative_error;
        result.analytic_sphere=true;
    }
    else
    {
        result.convexity_intervals=check_convexity(c);
        const Integrand integrand{c,std::hypot(direction[0],direction[2])/norm,std::abs(direction[1])/norm};
        std::priority_queue<Panel> panels;
        double value=0,error=0;
        for(unsigned i=0;i<32;++i)
        {
            const auto p=integrate_panel(integrand,pi*i/32,pi*(i+1)/32);
            panels.push(p);value+=p.value;error+=p.error;result.quadrature_evaluations+=15;
        }
        while(error>(std::max)(o.absolute_tolerance_dimensionless,o.relative_tolerance*std::abs(value)))
        {
            if(panels.size()>=o.maximum_quadrature_panels)
                throw std::runtime_error("Projected-area quadrature did not converge within its work limit.");
            const auto p=panels.top();panels.pop();const double mid=0.5*(p.lo+p.hi);
            if(mid==p.lo||mid==p.hi)throw std::runtime_error("Projected-area quadrature stagnated.");
            const auto l=integrate_panel(integrand,p.lo,mid),r=integrate_panel(integrand,mid,p.hi);
            value+=(l.value+r.value)-p.value;
            error=(std::max)(0.0,error-p.error)+l.error+r.error;
            panels.push(l);panels.push(r);result.quadrature_evaluations+=30;
        }
        // Re-sum final leaves; no reliance on long double being wider than double.
        double sum=0,correction=0,error_sum=0;
        while(!panels.empty())
        {
            const auto p=panels.top();panels.pop();
            const double y=p.value-correction,t=sum+y;correction=(t-sum)-y;sum=t;error_sum+=p.error;
        }
        if(error_sum>(std::max)(o.absolute_tolerance_dimensionless,o.relative_tolerance*std::abs(sum)))
            throw std::runtime_error("Projected-area final error check failed.");
        result.area_mm2=sum*radius2;
        result.estimated_absolute_error_mm2=error_sum*radius2;
        result.estimated_relative_error=error_sum/std::abs(sum);
    }
    if(!(result.area_mm2>0 && std::isfinite(result.area_mm2) &&
         std::isfinite(result.estimated_absolute_error_mm2)))
        throw std::runtime_error("Projected area is not representable.");
    return result;
}
}
