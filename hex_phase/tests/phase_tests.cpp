#include <ice_crystal/phase_generate.hpp>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace iceCrystal;
namespace
{
void check(bool x,const char* m){if(!x)throw std::runtime_error(m);}
void near(double x,double y,double eps,const char* m){check(finite_value(x)&&finite_value(y)&&std::fabs(x-y)<=eps,m);}
template<class F>void rejects(F&& f,const char* m){bool bad=false;try{f();}catch(const std::exception&){bad=true;}check(bad,m);}
std::vector<double> reconstruct(const HostCdf& c)
{
    std::vector<double> p(std::size_t(c.grid.nt)*c.grid.np);
    for(unsigned i=0;i<c.grid.nt;++i)for(unsigned j=0;j<c.grid.np;++j)
    {
        const auto k=std::size_t(j)*(c.grid.nt+1ull)+i;
        p[std::size_t(i)*c.grid.np+j]=(c.phi[j+1]-c.phi[j])*(c.theta[k+1]-c.theta[k]);
    }
    return p;
}
double l1(std::span<const double> a,std::span<const double> b)
{check(a.size()==b.size(),"array dimensions");Sum s;for(std::size_t i=0;i<a.size();++i)s.add(std::fabs(a[i]-b[i]));return s.value;}
void sha_plan()
{
    check(sha256("")=="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855","SHA256 empty");
    check(sha256("abc")=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad","SHA256 abc");
    check(sha256(std::string(1000000,'a'))=="cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0","SHA256 million a");
    const std::string header="qw,qx,qy,qz,weight,samples,seed\n";
    const auto p=parse_orientation_plan("# explicit\r\n"+header+"1,0,0,0,0.25,17,18446744073709551615\r\n1,0,0,0,0.75,55,3\n");
    check(p.total_samples==72&&p.nodes.size()==2,"orientation counts");near(p.nodes[1].weight,0.75,0,"number weight");
    rejects([&]{parse_orientation_plan(header+"1,0,0,0,2,1,1\n");},"unnormalized density accepted");
    rejects([&]{parse_orientation_plan(header+"2,0,0,0,1,1,1\n");},"nonunit quaternion accepted");
    rejects([&]{parse_orientation_plan(header+"1,0,0,0,nan,1,1\n");},"nan accepted");
    rejects([&]{parse_orientation_plan(header+"1,0,0,0,1,0,1\n");},"empty samples accepted");
    rejects([&]{parse_orientation_plan(header+"1,0,0,0,1,1,1,2\n");},"excess field accepted");
    rejects([&]{parse_orientation_plan(header+"1,0,0,0,0.5,18446744073709551615,1\n1,0,0,0,0.5,1,2\n");},"counter overflow accepted");
}
void rotations()
{
    const double v=std::sqrt(0.5);const auto r=rotation_from_quaternion(v,0,0,v);
    check((r.apply({1,0,0})-Vec3{0,1,0}).max_abs()<1e-15,"active rotation sign");
    check((r.inverse(r.apply({0.1,-3,7}))-Vec3{0.1,-3,7}).max_abs()<4e-15,"inverse");
    for(unsigned i=0;i<1000;++i)
    {
        const double angle=2*pi*sample_uniform(7,i,0);Vec3 a;
        check(normalize({2*sample_uniform(13,i,0)-1,2*sample_uniform(13,i,1)-1,2*sample_uniform(13,i,2)-1},a),"axis");
        const double s=std::sin(angle/2);const auto q=rotation_from_quaternion(std::cos(angle/2),s*a.x,s*a.y,s*a.z);
        near(q.apply({1,2,3}).dot(q.apply({1,2,3})),14,3e-14,"rotation norm");
    }
    Rotation bad;bad.z={0,0,-1};check(!bad.valid(),"reflection is not SO3");
}
void measure_contract()
{
    const auto g=make_phase_grid(64,128,{1,-0.4,0.3});std::vector<double> masses(64*128);
    for(unsigned i=0;i<64;++i)for(unsigned j=0;j<128;++j)masses[i*128+j]=(4*pi/128)*(g.edges[i+1]-g.edges[i]);
    const auto c=build_host_cdf(g,masses,64,128);near(c.report.g_stored,0,1e-15,"isotropic moment");
    for(unsigned j=0;j<=128;++j)near(c.phi[j],double(j)/128,1e-14,"uniform phi marginal");
    for(unsigned j=0;j<128;++j)for(unsigned i=0;i<=64;++i)near(c.theta[j*65+i],g.edges[i],2e-15,"uniform per-solid-angle conditional");
    const auto d=mass_to_density(g,masses);for(double v:d)near(v,1,3e-16,"divide mass by area");
    auto incorrect=masses;
    for(unsigned i=0;i<64;++i)for(unsigned j=0;j<128;++j)incorrect[i*128+j]*=g.edges[i+1]-g.edges[i];
    const auto wrong=reconstruct(build_host_cdf(g,incorrect,64,128));const auto right=reconstruct(c);
    check(l1(wrong,right)>0.2,"double area weighting test not sensitive");
    std::cout<<"  double_area_weighting_TV="<<0.5*l1(wrong,right)<<'\n';
}
void cells_and_poles()
{
    const auto g=make_phase_grid(8,16,{1,-0.4,0.3});PhaseMassAccumulator h(g);
    h.add(g.ki,3);h.add(-g.ki,1);const auto w=h.masses();
    for(unsigned j=0;j<16;++j){near(w[j],3./16,0,"forward cap");near(w[7*16+j],1./16,0,"backward cap");}
    const auto c=build_host_cdf(g,w,8,16);near(c.report.g_stored,0.5*(1-g.edges[1]),5e-16,"cell moment not ray moment");
    for(unsigned j=0;j<=16;++j)near(c.phi[j],double(j)/16,1e-15,"pole azimuth invariance");
    const auto r=rotation_from_quaternion(0.5,0.5,0.5,0.5);
    PhaseGrid rotated=g;rotated.ki=r.apply(g.ki);rotated.frame={r.apply(g.frame.e0),r.apply(g.frame.e1)};
    PhaseMassAccumulator other(rotated);
    for(unsigned i=0;i<8;++i)for(unsigned j=0;j<16;++j)
    {
        const double u=(g.edges[i]+g.edges[i+1])/2,mu=1-2*u,st=2*std::sqrt(u*(1-u)),phi=-pi+(j+.5)*2*pi/16;
        const auto d=g.ki*mu+g.frame.e0*(st*std::cos(phi))+g.frame.e1*(st*std::sin(phi));
        const auto b=phase_bin(g.view(),d);check(b.valid&&!b.pole&&b.index==i*16+j,"cell-center bin");
        h.add(d,0.01);other.add(r.apply(d),0.01);
    }
    other.add(rotated.ki,3);other.add(-rotated.ki,1);near(l1(h.masses(),other.masses()),0,2e-15,"frame covariance");
    const auto seam=phase_bin(g.view(),-g.frame.e0);check(seam.valid&&seam.index%16==0,"periodic seam");
    const auto axes=make_phase_grid(20,40,{0,0,1});
    for(double perturbation:{-1e-14,0.,1e-14})
    {
        const auto cell=phase_bin(axes.view(),{1,perturbation,0});
        check(cell.valid&&cell.index==10*40+20,"phi boundary roundoff assignment");
        const auto polar=phase_bin(axes.view(),{1,0,perturbation});
        check(polar.valid&&polar.index==10*40+20,"u boundary roundoff assignment");
    }
    Vec3 away;check(normalize({1,-1e-7,0.23},away),"boundary test vector");
    const auto unsnapped=phase_bin(axes.view(),away);
    check(unsnapped.index%40==19&&!unsnapped.boundary_snapped,"roundoff rule broadened into filter");
    rejects([&]{h.add({0,0,0},1);},"invalid direction accepted");
    rejects([&]{h.add(g.ki,0);},"zero weight accepted");
    rejects([&]{h.add(g.ki,std::numeric_limits<double>::infinity());},"infinite weight accepted");
}
void cdf_sparse_and_coarsen()
{
    const auto g=make_phase_grid(48,96,{0,-1,0});std::vector<double> mass(48*96);
    for(unsigned i=0;i<48;++i)for(unsigned j=0;j<96;++j)
        mass[i*96+j]=(j%7==0)?std::pow(10.,-12.*double((17*i+13*j)%97)/96):0;
    const auto c=build_host_cdf(g,mass,24,48);const auto p=reconstruct(c);
    near(c.report.cdf_mass,1,2e-15,"CDF normalized");check(c.report.l1_error<2e-14,"CDF roundoff");
    Sum total;for(double x:mass)total.add(x);std::vector<double> coarse(24*48);
    for(unsigned i=0;i<48;++i)for(unsigned j=0;j<96;++j)coarse[(i/2)*48+j/2]+=mass[i*96+j]/total.value;
    near(l1(coarse,p),0,2e-14,"conservative independent coarsening");
    const auto original=build_host_cdf(g,mass,48,96);
    for(unsigned j=0;j<96;++j)if(j%7)for(unsigned i=0;i<=48;++i)
        near(original.theta[j*49+i],g.edges[i],0,"zero marginal conditional fallback");
    CdfPolicy tight;tight.max_coarsening_tv=0;
    rejects([&]{build_host_cdf(g,mass,24,48,tight);},"coarsening quality ignored");
    rejects([&]{build_host_cdf(g,mass,25,48);},"noninteger coarsening accepted");
    std::vector<double> zero(mass.size());rejects([&]{build_host_cdf(g,zero,48,96);},"empty histogram accepted");
    mass[0]=-1;rejects([&]{build_host_cdf(g,mass,48,96);},"negative histogram accepted");
    const auto small=make_phase_grid(2,2,{1,0,0});std::vector<double> hdr{1,1e-20,1e-20,1e-20};
    CdfPolicy exact;exact.max_lost_mass=0;
    rejects([&]{build_host_cdf(small,hdr,2,2,exact);},"lost CDF mass silently repaired");
    std::cout<<"  sparse_CDF_L1="<<c.report.l1_error<<" coarsening_TV="<<c.report.coarsening_tv<<'\n';
}
struct Composed {std::vector<double> mass;CompositionAudit audit;};
Composed compose(PhaseGenerationSettings s,const OrientationPlan& plan,std::size_t chunk)
{
    PhaseMassAccumulator h(s.grid);CompositionAudit a;
    for(const auto& n:plan.nodes)
    {
        const double area=node_area(s,n);a.expected_area.add(n.weight*area);
        const auto c=node_trace_settings(s,n);const double scale=n.weight*area/n.samples;
        for(std::uint64_t first=0;first<n.samples;)
        {
            const auto count=std::min<std::uint64_t>(chunk,n.samples-first);const auto b=trace_batch_cpu(c,first,count);
            record_batch_audit(a,summarize_phase_rays(b.audits,s.max_balance_error),scale);h.add(b.outgoing,n.rotation,scale);first+=count;
        }
        ++a.orientations_done;
    }
    auto mass=h.masses();a.histogram=h.statistics();a.grid_moments=measure_mass(s.grid,mass);validate_composition(a,s,plan);
    return {mass,a};
}
void real_composition()
{
    const auto plan=parse_orientation_plan("qw,qx,qy,qz,weight,samples,seed\n1,0,0,0,0.25,73,12345\n"
        "0.7071067811865476,0,0,0.7071067811865476,0.75,211,54321\n");
    PhaseGenerationSettings s;s.prism={0.1,0.2};s.grid=make_phase_grid(36,72,{1,-0.4,0.3});
    s.out_nt=36;s.out_np=72;s.wavelength_nm=550;s.interior_index=1.31;
    validate_generation_settings(s);
    const auto a=compose(s,plan,10000),b=compose(s,plan,19);
    near(l1(a.mass,b.mass),0,1e-17,"batch partition changed histogram");
    const auto c=build_host_cdf(s.grid,a.mass,36,72);near(c.report.cdf_mass,1,1e-15,"composed CDF");
    // Each marginal is the SUM of unnormalized kernels, not a mean of already normalized CDFs.
    auto p0=plan,p1=plan;p0.nodes={plan.nodes[0]};p1.nodes={plan.nodes[1]};
    p0.nodes[0].weight=p1.nodes[0].weight=1;p0.total_samples=p0.nodes[0].samples;p1.total_samples=p1.nodes[0].samples;
    const auto one=compose(s,p0,43),two=compose(s,p1,43);std::vector<double> expected(a.mass.size());
    for(std::size_t i=0;i<expected.size();++i)expected[i]=0.25*one.mass[i]+0.75*two.mass[i];
    near(l1(expected,a.mass),0,3e-17,"orientation mixing weights");
    std::cout<<"  real_outputs="<<a.audit.trace.outputs<<" mass_error="<<a.audit.histogram_relative_error
             <<" unresolved="<<a.audit.unresolved_area.value/a.audit.expected_area.value<<" g="<<c.report.g_stored<<'\n';
    RayAudit invalid_interface{};invalid_interface.status=TraceStatus::Exhausted;
    invalid_interface.escaped_power=1;invalid_interface.unresolved_power=0;
    invalid_interface.max_interface_balance_error=std::numeric_limits<double>::quiet_NaN();
    check(!summarize_phase_rays({&invalid_interface,1}).accepted,"nonfinite interface diagnostic hidden");
    s.max_internal_hits=1;rejects([&]{compose(s,plan,100);},"truncation accepted");
    s.max_internal_hits=1024;s.tail_tolerance=1e-3;s.max_unresolved_fraction=1e-10;
    rejects([&]{compose(s,plan,100);},"global tail budget ignored");
}
void analytical_area_and_n1()
{
    const auto plan=parse_orientation_plan("qw,qx,qy,qz,weight,samples,seed\n1,0,0,0,0.5,17,11\n"
        "0.7071067811865476,0,0,0.7071067811865476,0.5,55,22\n");
    PhaseGenerationSettings s;s.prism={1,2};s.grid=make_phase_grid(16,32,{0,-1,0});s.out_nt=16;s.out_np=32;
    s.wavelength_nm=550;s.interior_index=1;s.tail_tolerance=0;
    const auto a=compose(s,plan,9);near(a.audit.expected_area.value,0.5*((3*sqrt3/2)+4),2e-15,"mean projected area");
    near(a.audit.escaped_area.value,a.audit.expected_area.value,4e-15,"transparent area");
    near(a.audit.histogram.forward_mass,a.audit.expected_area.value,4e-15,"forward mass retained");
    for(unsigned j=0;j<32;++j)near(a.mass[j],a.audit.expected_area.value/32,2e-16,"delta cap independent of basis");
    const auto c=build_host_cdf(s.grid,a.mass,16,32);near(c.report.g_stored,1-s.grid.edges[1],3e-15,"stored HG uses finite cells");
    // Independent numerical integration of the projection formula over S^2: E A = S/4.
    const HexPrism prism{0.2,0.7};Sum area;const unsigned n=100000;
    const double golden=pi*(3-std::sqrt(5.));
    for(unsigned i=0;i<n;++i)
    {
        const double mu=1-2*(i+0.5)/n,phi=golden*i,st=std::sqrt(1-mu*mu);
        area.add(prism.projected_area_mm2({st*std::cos(phi),mu,st*std::sin(phi)}));
    }
    const double relative=std::fabs(area.value/n/(prism.surface_area_mm2()/4)-1);
    check(relative<2e-6,"isotropic projected area S/4");std::cout<<"  area_S_over_4_relative="<<relative<<'\n';
}
}
int main()
{
    try
    {
        std::cout<<std::setprecision(17);
        const std::pair<const char*,std::function<void()>> tests[]{
            {"sha256_and_plan",sha_plan},{"rotations",rotations},{"probability_measure_contract",measure_contract},
            {"cells_poles_covariance",cells_and_poles},{"cdf_sparse_coarsening_rejection",cdf_sparse_and_coarsen},
            {"real_traces_composition_batch_invariance",real_composition},{"analytic_area_and_transparent_limit",analytical_area_and_n1}};
        for(const auto& [name,test]:tests){test();std::cout<<"PASS "<<name<<'\n';}
        return 0;
    }
    catch(const std::exception& e){std::cerr<<"FAILED: "<<e.what()<<'\n';return 1;}
}
