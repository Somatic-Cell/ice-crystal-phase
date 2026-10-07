#include <ice_crystal/hex_trace.hpp>
#include <algorithm>
#include <array>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace iceCrystal;
namespace
{
void check(bool b,const char* m){if(!b)throw std::runtime_error(m);}
void near(double a,double b,double t,const char* message)
{
    if(!finite_value(a)||!finite_value(b)||::fabs(a-b)>t)
        throw std::runtime_error(std::string(message)+": "+std::to_string(a)+" vs "+std::to_string(b));
}
Vec3 unit(Vec3 v){Vec3 r{};check(normalize(v,r),"normalize");return r;}
void near_vec(Vec3 a,Vec3 b,double t,const char* m){check((a-b).max_abs()<=t,m);}
TraceSettings config(Vec3 k={0,-1,0},std::uint64_t n=1024,unsigned hits=256,double tol=1e-12)
{return make_trace_settings({1,2},k,1,1.31,n,12345,hits,tol);}
template<class F> void rejects(F f)
{bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}check(rejected,"input not rejected");}
Vec3 rotate60(Vec3 p){return {0.5*p.x-sqrt3/2*p.z,p.y,sqrt3/2*p.x+0.5*p.z};}

void projection()
{
    const HexPrism h{1,2};
    near(h.projected_area_mm2({0,1,0}),3*sqrt3/2,1e-14,"axial area");
    near(h.projected_area_mm2({1,0,0}),4,1e-14,"side area");
    Sum areas{};
    constexpr unsigned n=100000;
    for(unsigned i=0;i<n;++i)
    {
        const double y=1-2*(i+0.5)/n,phi=2*pi*i*0.6180339887498948482;
        const double r=::sqrt((1-y)*(1+y));Vec3 k{r*::cos(phi),y,r*::sin(phi)};
        const double a=h.projected_area_mm2(k);
        near(a,h.projected_area_mm2(-k),5e-14,"central symmetry of projected area");
        near(a,h.projected_area_mm2(rotate60(k)),5e-14,"60 degree symmetry of area");
        areas.add(a);
    }
    near(areas.value/n,h.surface_area_mm2()/4,2e-5,"Cauchy mean projected area");
}
// Independent Moller-Trumbore mesh oracle. Only test code uses triangles;
// production prism intersections do not use a mesh or BVH.
bool triangle_hit(Vec3 p,Vec3 d,Vec3 a,Vec3 b,Vec3 c,double& t)
{
    Vec3 e1=b-a,e2=c-a,h=d.cross(e2);double det=e1.dot(h);
    if(::fabs(det)<1e-13)return false;
    Vec3 s=p-a;double u=s.dot(h)/det;if(u<0||u>1)return false;
    Vec3 q=s.cross(e1);double v=d.dot(q)/det;if(v<0||u+v>1)return false;
    t=e2.dot(q)/det;return t>=0;
}
void intersection_mesh()
{
    const HexPrism h{1,2};
    std::vector<std::array<Vec3,3>> triangles;
    for(unsigned f=0;f<6;++f)
    {
        Vec3 a=h.vertex(f),b=h.vertex(f+1);a.y=1;b.y=1;
        Vec3 c=a,d=b;c.y=-1;d.y=-1;
        triangles.push_back({a,b,c});triangles.push_back({b,d,c});
        triangles.push_back({Vec3{0,1,0},a,b});triangles.push_back({Vec3{0,-1,0},c,d});
    }
    double largest=0;
    for(unsigned i=0;i<10000;++i)
    {
        Vec3 p{6*sample_uniform(47,i,0)-3,6*sample_uniform(47,i,1)-3,6*sample_uniform(47,i,2)-3};
        Vec3 d=unit(Vec3{sample_uniform(49,i,0)-0.5,sample_uniform(49,i,1)-0.5,sample_uniform(49,i,2)-0.5});
        double best=DBL_MAX;
        for(const auto& v:triangles){double t=0;if(triangle_hit(p,d,v[0],v[1],v[2],t))best=(std::min)(best,t);}
        const auto hit=intersect_prism(h,p,d);
        if(best==DBL_MAX){check(hit.status==HitStatus::Miss,"mesh/clipping miss disagreement");continue;}
        check(hit.status==HitStatus::Regular,"unexpected ambiguous random hit");
        near(hit.distance,best,2e-12,"mesh/clipping distance");
        largest=(std::max)(largest,::fabs(hit.distance-best));
    }
    std::cout<<"  maximum_mesh_distance_error="<<largest<<'\n';
    auto edge=intersect_prism(h,{3,0,0},unit({-1,0,0}));
    check(edge.status==HitStatus::Regular,"face hit");
    edge=intersect_prism(h,{0,0,3},{0,0,-1});
    check(edge.status==HitStatus::AmbiguousBoundary,"vertex normal must be ambiguous");
}
void entry_sampling()
{
    auto c=config(unit({1,-0.7,0.22}));const auto& h=c.prism;const auto k=c.incident_direction;
    constexpr unsigned n=100000;std::uint64_t counts[8]{};Sum radial{};
    for(unsigned i=0;i<n;++i)
    {
        const double u0=sample_uniform(55,i,0),u1=sample_uniform(55,i,1),u2=sample_uniform(55,i,2);
        const auto s=sample_entry(h,k,u0,u1,u2);
        check(s.valid&&strict_face_point(h,s.position,s.face),"sample outside strict face");++counts[s.face];
        const auto hit=intersect_prism(h,s.position-k*6,k);
        check(hit.status==HitStatus::Regular&&hit.face==s.face,"first hit not sampled face");
        near_vec(hit.position,s.position,2e-12,"first-hit point mismatch");
        const auto cap=sample_entry(h,{0,-1,0},u0,u1,u2);
        check(cap.face==6,"cap selector");radial.add(cap.position.x*cap.position.x+cap.position.z*cap.position.z);
    }
    const double a=h.projected_area_mm2(k);
    for(unsigned f=0;f<8;++f)
    {
        const double p=h.face_area_mm2(f)*max_value(0,-h.normal(f).dot(k))/a;
        near(double(counts[f]),n*p,7*::sqrt(n*p*(1-p))+3,"projected face frequencies");
    }
    near(radial.value/n,5.0/12,0.004,"uniform hexagonal cap second moment");
    for(unsigned i=0;i<10000;++i)
        for(unsigned j=0;j<3;++j)check(sample_uniform(0,i,j)>0&&sample_uniform(0,i,j)<1,"open RNG");
}
FluxRay plane_ray(double angle)
{
    FluxRay r{};r.direction={::sin(angle),-::cos(angle),0};r.frame={{0,0,1},r.direction.cross({0,0,1})};
    r.response=JonesFlux::identity();return r;
}
void fresnel()
{
    auto r=plane_ray(0);const auto normal=Vec3{0,1,0};
    const double R=(1-1.31)*(1-1.31)/((1+1.31)*(1+1.31));
    auto b=split_flux(r,normal,1,1.31);check(b.valid&&b.has_transmission,"normal Fresnel");
    near(b.reflected.response.power(),R,2e-15,"normal R");
    near(b.transmitted.response.power(),1-R,2e-15,"normal T power factor");
    double worst=0;
    for(unsigned i=0;i<10000;++i)
    {
        const double angle=(i+0.5)/10000*pi/2;
        for(unsigned reverse=0;reverse<2;++reverse)
        {
            const double ni=reverse?1.31:1,nt=reverse?1:1.31;
            auto incoming=plane_ray(angle);auto out=split_flux(incoming,normal,ni,nt);
            check(out.valid,"interface rejected in angle sweep");
            const double total=out.reflected.response.power()+(out.has_transmission?out.transmitted.response.power():0);
            near(total,1,2e-13,"interface energy conservation");
            worst=(std::max)(worst,::fabs(total-1));
            if(out.has_transmission)
                near(ni*::sin(angle),nt*out.transmitted.direction.x,3e-14,"Snell tangent");
        }
    }
    auto brewster=plane_ray(::atan(1.31));brewster.response={};brewster.response.column[0].y.re=::sqrt(2.0);
    b=split_flux(brewster,normal,1,1.31);
    check(b.valid,"Brewster interface");check(b.reflected.response.power()<1e-28,"Brewster p reflection");
    auto tir=plane_ray(70*pi/180);b=split_flux(tir,normal,1.31,1);
    check(b.valid&&b.total_internal_reflection&&!b.has_transmission,"TIR branch");
    near(b.reflected.response.power(),1,2e-14,"TIR norm");
    check(::fabs(b.reflected.response.column[0].x.im)>0.1,"TIR complex phase retained");
    const auto phase_s=b.reflected.response.column[0].x,phase_p=b.reflected.response.column[1].y;
    check(::fabs(phase_s.re*phase_p.im-phase_s.im*phase_p.re)>0.01,"TIR relative retardance retained");
    std::cout<<"  maximum_interface_power_error="<<worst<<'\n';
}
void polarization_history()
{
    auto in=plane_ray(55*pi/180);
    auto b1=split_flux(in,{0,1,0},1,1.31);
    auto b2=split_flux(b1.transmitted,{0,1,0},1.31,1);
    check(b1.valid&&b2.valid,"two-surface chain");
    const double ts=b1.transmitted.response.column[0].x.norm2();
    const double tp=b1.transmitted.response.column[1].y.norm2();
    near(b2.transmitted.response.power(),0.5*(ts*ts+tp*tp),1e-14,"polarization correlation through slab");
    const double reset_every_surface=0.25*(ts+tp)*(ts+tp);
    check(::fabs(b2.transmitted.response.power()-reset_every_surface)>1e-4,"not reset to unpolarized at every interface");
    const double angle=0.43;
    Frame turned{in.frame.e0*::cos(angle)+in.frame.e1*::sin(angle),
                 in.frame.e1*::cos(angle)-in.frame.e0*::sin(angle)};
    auto other=in;other.response=change_frame(in.response,in.frame,turned);other.frame=turned;
    auto b3=split_flux(other,{0,1,0},1,1.31);
    near(b3.transmitted.response.power(),b1.transmitted.response.power(),2e-14,"basis covariance");
}
void slab_and_tail()
{
    auto c=config({0,-1,0},10,2,1e-15);
    const double R=(1-1.31)*(1-1.31)/((1+1.31)*(1+1.31));
    auto batch=trace_batch_cpu(c,0,10);const auto s=summarize(batch.audits);
    check(!s.accepted&&s.status_counts[2]==10,"truncated rays must not pass");
    for(unsigned i=0;i<10;++i)
    {
        const auto& a=batch.audits[i];check(a.output_count==3,"R + two transmitted exits");
        near(a.escaped_power,R+(1-R)*(1-R)*(1+R),3e-14,"analytic slab escaped");
        near(a.unresolved_power,(1-R)*R*R,3e-15,"analytic slab residual");
        near(a.balance_error,0,3e-14,"slab total");
        const auto& direct=batch.outgoing[batch.offsets[i]+1];
        near(direct.internal_length_mm,2,1e-14,"direct slab length");
        near_vec(direct.direction,{0,-1,0},2e-14,"forward delta retained");
    }
    rejects([&]{require_accepted(s);});
    c.max_internal_hits=256;c.residual_power_tolerance=1e-12;
    batch=trace_batch_cpu(c,0,10);require_accepted(summarize(batch.audits));
    for(const auto& a:batch.audits)check(a.unresolved_power<=c.residual_power_tolerance,"tail bound");
    c.interior_index=1;c.residual_power_tolerance=0;
    batch=trace_batch_cpu(c,0,10);require_accepted(summarize(batch.audits));
    for(const auto& a:batch.audits){check(a.output_count==1,"index-matched event count");near(a.escaped_power,1,1e-15,"index matched power");}
}
void prism_trace_and_batch()
{
    auto c=config(unit({1,-0.4,0.3}),2048);
    auto batch=trace_batch_cpu(c,0,2048);auto sum=summarize(batch.audits);require_accepted(sum);
    double max_error=0;
    for(std::size_t i=0;i<batch.audits.size();++i)
    {
        Sum power{};
        for(auto j=batch.offsets[i];j<batch.offsets[i+1];++j)
        {
            const auto& o=batch.outgoing[j];power.add(o.power_fraction);
            check(o.frame.valid(o.direction),"output frame");
            check(o.exit_face<8&&o.entry_face<8,"output faces");
            check(c.prism.normal(o.exit_face).dot(o.direction)>=-1e-14,"outgoing ray points inward");
            near(o.power_fraction,o.response.power(),2e-14,"Jones/weight consistency");
            check(strict_face_point(c.prism,o.position_mm/c.prism.circumradius_mm,o.exit_face),"exit on facet");
        }
        near(power.value,batch.audits[i].escaped_power,1e-14,"event sum audit");
        max_error=(std::max)(max_error,::fabs(batch.audits[i].balance_error));
    }
    auto a=trace_batch_cpu(c,0,333),b=trace_batch_cpu(c,333,2048-333);
    check(a.outgoing.size()+b.outgoing.size()==batch.outgoing.size(),"batch-invariant output count");
    for(std::size_t i=0;i<b.outgoing.size();++i)
    {
        const auto& u=b.outgoing[i];const auto& v=batch.outgoing[a.outgoing.size()+i];
        check(u.power_fraction==v.power_fraction&&u.incident_sample_id==v.incident_sample_id,"batch-invariant sample");
        near_vec(u.direction,v.direction,0,"batch-invariant direction");
    }
    std::cout<<"  oblique_outputs="<<sum.outputs<<" maximum_ray_balance_error="<<max_error
        <<" mean_unresolved="<<sum.unresolved_power_sum/sum.rays<<'\n';
}
void symmetry_and_scale()
{
    auto c=config(unit({1,-0.4,0.3}),100);auto d=c;d.incident_direction=rotate60(c.incident_direction);
    check(make_frame(d.incident_direction,d.incident_frame),"rotated frame");
    for(unsigned i=0;i<100;++i)
    {
        auto entry=sample_entry(c.prism,c.incident_direction,sample_uniform(43,i,0),sample_uniform(43,i,1),sample_uniform(43,i,2));
        auto rotated=entry;rotated.position=rotate60(entry.position);
        if(rotated.face<6)rotated.face=(rotated.face+1)%6;
        std::array<OutgoingSample,300> oa{},ob{};BufferSink sa{oa.data(),oa.size()},sb{ob.data(),ob.size()};
        const auto a=trace_from_entry(c,i,entry,sa),b=trace_from_entry(d,i,rotated,sb);
        check(accepted_status(a.status)&&accepted_status(b.status),"rotational trace accepted");
        check(sa.size==sb.size,"rotational event count");
        for(std::size_t j=0;j<sa.size;++j)
        { near(oa[j].power_fraction,ob[j].power_fraction,2e-12,"rotational power");near_vec(rotate60(oa[j].direction),ob[j].direction,2e-12,"rotational direction"); }
    }
    const auto base=trace_batch_cpu(c,0,100);
    for(double scale:{1e-6,1e6})
    {
        auto s=c;s.prism.circumradius_mm*=scale;s.prism.length_mm*=scale;
        auto b=trace_batch_cpu(s,0,100);require_accepted(summarize(b.audits));
        check(base.outgoing.size()==b.outgoing.size(),"scale event count");
        for(std::size_t j=0;j<b.outgoing.size();++j)
        {
            near(base.outgoing[j].power_fraction,b.outgoing[j].power_fraction,1e-12,"scale power");
            near_vec(base.outgoing[j].direction,b.outgoing[j].direction,1e-12,"scale direction");
            near(base.outgoing[j].internal_length_mm,b.outgoing[j].internal_length_mm/scale,1e-10,"scale path length");
        }
    }
}
void direction_aspect_sweep()
{
    std::uint64_t statuses[8]{},tir=0;
    double max_balance=0,largest_tail=0;
    for(unsigned condition=0;condition<300;++condition)
    {
        const Vec3 k{2*sample_uniform(119,condition,0)-1,
                     2*sample_uniform(119,condition,1)-1,
                     2*sample_uniform(119,condition,2)-1};
        const double length=0.2+5*sample_uniform(127,condition,0);
        auto c=make_trace_settings({1,length},k,1,1.31,100,condition+700,1024,1e-12);
        for(unsigned i=0;i<100;++i)
        {
            CountSink sink{};const auto a=trace_sample(c,i,sink);
            ++statuses[static_cast<unsigned>(a.status)];tir+=a.tir_count;
            max_balance=(std::max)(max_balance,::fabs(a.balance_error));
            largest_tail=(std::max)(largest_tail,a.unresolved_power);
            if(!accepted_status(a.status))
                std::cerr<<"sweep failure condition="<<condition<<" sample="<<i<<" status="<<status_name(a.status)
                    <<" hits="<<a.internal_hits<<" escaped="<<a.escaped_power<<" residual="<<a.unresolved_power
                    <<" length="<<length<<" k="<<k.x<<","<<k.y<<","<<k.z<<'\n';
            check(accepted_status(a.status),"direction/aspect sweep unaccepted ray");
            check(::fabs(a.balance_error)<1e-10,"direction/aspect sweep balance");
        }
    }
    check(tir>0,"sweep did not exercise total internal reflection");
    std::cout<<"  sweep_rays=30000 tir_encounters="<<tir<<" maximum_balance_error="<<max_balance
        <<" maximum_unresolved="<<largest_tail<<'\n';
}
void long_tir_residual()
{
    const unsigned condition=157;
    const Vec3 k{2*sample_uniform(119,condition,0)-1,
                 2*sample_uniform(119,condition,1)-1,
                 2*sample_uniform(119,condition,2)-1};
    const double length=0.2+5*sample_uniform(127,condition,0);
    auto c=make_trace_settings({1,length},k,1,1.31,100,condition+700,256,1e-12);
    CountSink sink{};const auto limited=trace_sample(c,93,sink);
    check(limited.status==TraceStatus::InteractionLimit,"long TIR must report cap");
    check(limited.unresolved_power>0.98&&limited.tir_count==256,"large missing flux must remain visible");
    c.max_internal_hits=1024;
    const auto resolved=trace_sample(c,93,sink);
    check(accepted_status(resolved.status)&&resolved.internal_hits>256,"higher-order trace resolves escaped ray");
    check(resolved.unresolved_power<1e-12,"long TIR final residual");
    near(resolved.escaped_power+resolved.unresolved_power,1,1e-10,"long TIR total balance");
    std::cout<<"  long_TIR_hits="<<resolved.internal_hits<<" capped_residual="<<limited.unresolved_power
        <<" resolved_residual="<<resolved.unresolved_power<<'\n';
}
void errors()
{
    rejects([]{auto s=config();s.prism.circumradius_mm=0;validate_settings(s);});
    rejects([]{auto s=config();s.interior_index=0;validate_settings(s);});
    rejects([]{auto s=config();s.incident_direction={0,0,0};validate_settings(s);});
    rejects([]{auto s=config();s.residual_power_tolerance=1;validate_settings(s);});
    rejects([]{auto s=config();s.max_internal_hits=0;validate_settings(s);});
    rejects([]{auto s=config();s.total_incident_samples=0;validate_settings(s);});
    rejects([]{static_cast<void>(trace_batch_cpu(config(),1024,1));});
    rejects([]{static_cast<void>(trace_batch_cpu(config(),0,1,0));});
    auto c=config();EntrySample vertex{{0,1,1},6,true};CountSink counter{};
    const auto a=trace_from_entry(c,0,vertex,counter);
    check(a.status==TraceStatus::AmbiguousBoundary&&a.unresolved_power==1,"ambiguous source not discarded silently");
    auto e=sample_entry(c.prism,c.incident_direction,0.3,0.4,0.6);BufferSink sink{};
    const auto b=trace_from_entry(c,0,e,sink);
    check(b.status==TraceStatus::OutputOverflow&&b.unresolved_power==1,"bounded output");
    check(!sample_entry(c.prism,c.incident_direction,1,0.3,0.3).valid,"endpoint sampling rejected");
    std::array<RayAudit,1> invalid{};invalid[0].status=TraceStatus::Exhausted;
    invalid[0].escaped_power=std::numeric_limits<double>::quiet_NaN();
    check(!summarize(invalid).accepted,"nonfinite output rejected");
}
}
int main()
{
    const std::array<std::pair<const char*,std::function<void()>>,11> tests{{
        {"projected_area",projection},{"independent_mesh_intersection",intersection_mesh},
        {"projected_entry_sampling",entry_sampling},{"fresnel_snell_tir",fresnel},
        {"polarization_history",polarization_history},{"analytic_slab_and_tail",slab_and_tail},
        {"oblique_trace_and_batching",prism_trace_and_batch},{"symmetry_and_scale",symmetry_and_scale},
        {"direction_aspect_sweep",direction_aspect_sweep},{"long_TIR_residual",long_tir_residual},{"failure_paths",errors}}};
    unsigned passed=0;std::cout<<std::setprecision(17);
    for(const auto& [name,test]:tests)
    {
        try{test();++passed;std::cout<<"PASS "<<name<<'\n';}
        catch(const std::exception& e){std::cerr<<"FAIL "<<name<<": "<<e.what()<<'\n';return 1;}
    }
    std::cout<<passed<<"/"<<tests.size()<<" test groups passed\n";return 0;
}
