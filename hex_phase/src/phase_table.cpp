#include <ice_crystal/phase_table.hpp>
#include <rainbow/phase_cdf_math.hpp>
#include <rainbow/npy_writer.hpp>
#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>

namespace iceCrystal
{
Rotation rotation_from_quaternion(double w,double x,double y,double z)
{
    const double n=w*w+x*x+y*y+z*z;
    if(!finite_value(n)||::fabs(n-1)>2e-12)throw std::invalid_argument("Quaternion must be unit length (w,x,y,z).");
    const double k=1/std::sqrt(n);w*=k;x*=k;y*=k;z*=k;
    Rotation r{{1-2*(y*y+z*z),2*(x*y+w*z),2*(x*z-w*y)},
               {2*(x*y-w*z),1-2*(x*x+z*z),2*(y*z+w*x)},
               {2*(x*z+w*y),2*(y*z-w*x),1-2*(x*x+y*y)}};
    if(!r.valid())throw std::invalid_argument("Quaternion rotation failed orthonormality check.");
    return r;
}
PhaseGrid make_phase_grid(std::uint32_t nt,std::uint32_t np,Vec3 ki)
{
    if(!nt||!np||std::uint64_t(nt)*np>0xffffffffull||nt>32u*65535u||np>0x7fffffffu)
        throw std::invalid_argument("Unsupported phase grid dimensions.");
    PhaseGrid g;g.nt=nt;g.np=np;
    if(!normalize(ki,g.ki)||!make_frame(g.ki,g.frame))throw std::invalid_argument("Invalid incident direction.");
    g.edges=rainbow::make_phase_u_edges(nt);validate_phase_grid(g);return g;
}
void validate_phase_grid(const PhaseGrid& g)
{
    if(!g.nt||!g.np||g.edges.size()!=std::size_t(g.nt)+1||!g.frame.valid(g.ki)||
       g.edges.front()!=0||g.edges.back()!=1||std::uint64_t(g.nt)*g.np>0xffffffffull)
        throw std::invalid_argument("Malformed phase grid.");
    for(std::uint32_t i=0;i<g.nt;++i)
        if(!finite_value(g.edges[i+1])||!(g.edges[i+1]>g.edges[i]))
            throw std::invalid_argument("Nonincreasing u edges.");
}
PhaseMassAccumulator::PhaseMassAccumulator(PhaseGrid g):grid_(std::move(g))
{validate_phase_grid(grid_);mass_.resize(std::size_t(grid_.nt)*grid_.np);}
void PhaseMassAccumulator::add(Vec3 direction,double w)
{
    if(!(w>0)||!finite_value(w))throw std::runtime_error("Histogram requires positive finite weight; no lost sample is repaired.");
    const auto b=phase_bin(grid_.view(),direction);
    if(!b.valid)throw std::runtime_error("Invalid outgoing direction during accumulation.");
    if(count_==(std::numeric_limits<std::uint64_t>::max)())throw std::overflow_error("Output counter overflow.");
    ++count_;snapped_+=b.boundary_snapped;point_mass_.add(w);point_axial_.add(w*b.mu);
    if(b.pole>0)forward_.add(w);
    else if(b.pole<0)backward_.add(w);
    else mass_[b.index].add(w);
}
void PhaseMassAccumulator::add(std::span<const OutgoingSample> samples,const Rotation& r,double scale)
{
    if(!r.valid()||!(scale>0)||!finite_value(scale))throw std::invalid_argument("Invalid orientation / area weight.");
    for(const auto& s:samples)
    {
        const double w=scale*s.power_fraction;
        if(!(s.power_fraction>0)||!finite_value(s.power_fraction)||!(w>0)||!finite_value(w))
            throw std::runtime_error("Output area weight over/underflow.");
        add(r.apply(s.direction),w);
    }
}
std::vector<double> PhaseMassAccumulator::masses() const
{
    std::vector<double> result(mass_.size());
    const double f=forward_.value/grid_.np,b=backward_.value/grid_.np;
    if((forward_.value>0&&!(f>0))||(backward_.value>0&&!(b>0)))throw std::runtime_error("Pole weight underflow.");
    for(std::size_t k=0;k<result.size();++k)
    {
        Sum s=mass_[k];const auto i=k/grid_.np;
        if(i==0)s.add(f);
        if(i+1==grid_.nt)s.add(b);
        result[k]=s.value;
    }
    return result;
}
HistogramStats PhaseMassAccumulator::statistics() const
{return {point_mass_.value,point_axial_.value,forward_.value,backward_.value,0,count_,0,0,snapped_};}
HistogramMoments measure_mass(const PhaseGrid& g,std::span<const double> mass)
{
    validate_phase_grid(g);
    if(mass.size()!=std::uint64_t(g.nt)*g.np)throw std::invalid_argument("Mass array dimensions differ.");
    Sum total{},axial{};
    for(std::size_t k=0;k<mass.size();++k)
    {
        const double w=mass[k];
        if(!finite_value(w)||w<0)throw std::invalid_argument("Negative/nonfinite mass.");
        const auto i=k/g.np;total.add(w);axial.add(w*(1-g.edges[i]-g.edges[i+1]));
    }
    if(!(total.value>0)||!finite_value(total.value)||!finite_value(axial.value))
        throw std::runtime_error("Empty/nonfinite histogram.");
    return {total.value,axial.value};
}
std::vector<double> mass_to_density(const PhaseGrid& g,std::span<const double> mass)
{
    static_cast<void>(measure_mass(g,mass));std::vector<double> d(mass.size());
    for(std::size_t k=0;k<mass.size();++k)
    {
        const auto i=k/g.np;d[k]=mass[k]/((4*pi/g.np)*(g.edges[i+1]-g.edges[i]));
        if(!finite_value(d[k])||(mass[k]>0&&!(d[k]>0)))throw std::runtime_error("Mass to density over/underflow.");
    }
    return d;
}
void validate_cdf_policy(const CdfPolicy& p)
{
    if(!finite_value(p.max_l1_error)||p.max_l1_error<0||p.max_l1_error>2||
       !finite_value(p.max_lost_mass)||p.max_lost_mass<0||p.max_lost_mass>1||
       !finite_value(p.max_coarsening_tv)||p.max_coarsening_tv<0||p.max_coarsening_tv>1)
        throw std::invalid_argument("Invalid CDF quality policy.");
}
HostCdf build_host_cdf(const PhaseGrid& g,std::span<const double> mass,
    std::uint32_t nt,std::uint32_t np,const CdfPolicy& policy)
{
    validate_cdf_policy(policy);
    const auto fine=measure_mass(g,mass);
    if(!nt||!np||nt>g.nt||np>g.np||g.nt%nt||g.np%np)
        throw std::invalid_argument("Storage dimensions must divide histogram dimensions.");
    HostCdf out;out.grid=make_phase_grid(nt,np,g.ki);out.grid.frame=g.frame;
    std::vector<Sum> sums(std::size_t(nt)*np);
    const auto dt=g.nt/nt,dp=g.np/np;
    for(std::uint32_t i=0;i<=nt;++i)out.grid.edges[i]=g.edges[std::size_t(i)*dt];
    for(std::uint32_t i=0;i<g.nt;++i)for(std::uint32_t j=0;j<g.np;++j)
        sums[std::size_t(i/dt)*np+j/dp].add(mass[std::size_t(i)*g.np+j]);
    std::vector<double> coarse(sums.size());
    for(std::size_t k=0;k<coarse.size();++k)coarse[k]=sums[k].value;
    const auto stored=measure_mass(out.grid,coarse);
    auto& report=out.report;report.g_source=fine.axial/fine.mass;
    report.aggregation_relative_error=std::fabs(stored.mass/fine.mass-1);
    Sum tv{};
    for(std::uint32_t i=0;i<g.nt;++i)for(std::uint32_t j=0;j<g.np;++j)
    {
        const auto ci=i/dt,cj=j/dp;
        const double fraction=(g.edges[i+1]-g.edges[i])/
            (out.grid.edges[ci+1]-out.grid.edges[ci])/dp;
        tv.add(std::fabs(mass[std::size_t(i)*g.np+j]/fine.mass-
            coarse[std::size_t(ci)*np+cj]/stored.mass*fraction));
    }
    report.coarsening_tv=0.5*tv.value;
    if(!finite_value(report.coarsening_tv)||report.coarsening_tv>policy.max_coarsening_tv||
       report.aggregation_relative_error>1e-12)throw std::runtime_error("Conservative coarsening audit failed.");
    // Ordered, nonnegative FP64 scans retain plateaus. Do not repair a CDF.
    const double largest=*std::max_element(coarse.begin(),coarse.end());
    out.phi.assign(std::size_t(np)+1,0);out.theta.assign(std::size_t(np)*(nt+1ull),0);
    std::vector<double> columns(np);
    for(std::uint32_t j=0;j<np;++j)
    {
        const auto base=std::size_t(j)*(nt+1ull);double s=0;
        for(std::uint32_t i=0;i<nt;++i)
        {
            const double w=coarse[std::size_t(i)*np+j]/largest;
            if(coarse[std::size_t(i)*np+j]>0&&!(w>0))throw std::runtime_error("CDF weight underflow.");
            s+=w;out.theta[base+i+1]=s;
        }
        columns[j]=s;out.phi[j+1]=out.phi[j]+s;
    }
    const double total=out.phi.back();
    if(!(total>0)||!finite_value(total))throw std::runtime_error("Invalid CDF total.");
    for(std::uint32_t j=0;j<np;++j)
    {
        const auto base=std::size_t(j)*(nt+1ull);
        for(std::uint32_t i=1;i<nt;++i)
            out.theta[base+i]=columns[j]>0?out.theta[base+i]/columns[j]:out.grid.edges[i];
        out.theta[base+nt]=1;
    }
    for(std::uint32_t j=1;j<np;++j)out.phi[j]/=total;
    out.phi.back()=1;
    Sum l1{},lost{},recovered{},axial{};
    for(std::uint32_t j=0;j<np;++j)for(std::uint32_t i=0;i<nt;++i)
    {
        const auto base=std::size_t(j)*(nt+1ull)+i;
        const double q=(out.phi[j+1]-out.phi[j])*(out.theta[base+1]-out.theta[base]);
        const double p=coarse[std::size_t(i)*np+j]/stored.mass;
        if(!finite_value(q)||q<0||out.phi[j+1]<out.phi[j]||out.theta[base+1]<out.theta[base])
            throw std::runtime_error("Malformed CDF; refusing repair.");
        if(coarse[std::size_t(i)*np+j]>0&&!(p>0))throw std::runtime_error("Normalized CDF weight underflow.");
        const double error=std::fabs(p-q);l1.add(error);
        report.maximum_cell_error=std::max(report.maximum_cell_error,error);
        if(p>0&&q==0){lost.add(p);++report.lost_cells;}
        recovered.add(q);axial.add(q*(1-out.grid.edges[i]-out.grid.edges[i+1]));
    }
    report.l1_error=l1.value;report.lost_mass=lost.value;report.cdf_mass=recovered.value;
    report.g_stored=axial.value/recovered.value;
    if(report.l1_error>policy.max_l1_error||report.lost_mass>policy.max_lost_mass||
       std::fabs(report.cdf_mass-1)>1e-10||!finite_value(report.g_stored)||
       !(report.g_stored>-1&&report.g_stored<1))throw std::runtime_error("CDF mass / HG audit failed.");
    return out;
}
void HostCdf::write_arrays(const std::filesystem::path& p) const
{
    const std::array<std::uint64_t,1> ashape{std::uint64_t(grid.np)+1},ushape{std::uint64_t(grid.nt)+1};
    const std::array<std::uint64_t,2> tshape{grid.np,std::uint64_t(grid.nt)+1};
    rainbow::NpyFloat64Writer a(p/"phi_cdf.npy",ashape);a.append(phi);a.finish();
    rainbow::NpyFloat64Writer t(p/"theta_given_phi_cdf.npy",tshape);t.append(theta);t.finish();
    rainbow::NpyFloat64Writer u(p/"u_edges.npy",ushape);u.append(grid.edges);u.finish();
}
} // namespace iceCrystal
