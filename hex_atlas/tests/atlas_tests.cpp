#include <ice_crystal/atlas.hpp>
#include <algorithm>
#include <array>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace iceCrystal;
namespace
{
void require(bool b,const char* why){if(!b)throw std::runtime_error(why);}
template<class F> void rejects(F f){bool hit=false;try{f();}catch(const std::exception&){hit=true;}require(hit,"Expected rejection");}
atlas::Settings settings()
{
    atlas::Settings s;s.prism={.1,.2};s.interior_index=1.31;s.wavelength_nm=550;
    s.incident_theta_intervals=3;s.incident_phi_intervals=2;s.theta_count=12;s.phi_count=24;s.area_subdivisions=2;return s;
}
Vec3 transform(Vec3 v,unsigned rot,bool mirror,bool flip)
{
    if(mirror)v.z=-v.z;
    if(flip)v.y=-v.y;
    const double a=rot*pi/3,c=std::cos(a),s=std::sin(a);
    return {c*v.x-s*v.z,v.y,s*v.x+c*v.z};
}
struct VectorSink
{
    std::vector<OutgoingSample> samples;
    bool push(const OutgoingSample& s){samples.push_back(s);return true;}
};
void area_test()
{
    const HexPrism h{.1,.7};
    for(Vec3 k: {Vec3{1,0,0},Vec3{0,1,0},Vec3{1,2,3},Vec3{-1,-.001,-.9}})
    {
        normalize(k,k);const auto r=atlas::make_area_rule(h,k,3);Sum sum{};
        for(const auto& a:r){require(a.projected_weight_mm2>0,"Nonpositive area");require(strict_face_point(h,a.entry.position,a.entry.face),"Invalid area point");sum.add(a.projected_weight_mm2);}
        require(std::fabs(sum.value/h.projected_area_mm2(k)-1)<5e-14,"Area integration");
        const auto same=atlas::make_area_rule(h,k,3);require(same.size()==r.size(),"Nondeterministic count");
        for(std::size_t i=0;i<r.size();++i)require((same[i].entry.position-r[i].entry.position).max_abs()==0 && same[i].projected_weight_mm2==r[i].projected_weight_mm2,"Nondeterministic quadrature");
    }
}
void nodes_test()
{
    auto s=settings();
    for(unsigned i=0;i<=s.incident_theta_intervals;++i)for(unsigned j=0;j<=s.incident_phi_intervals;++j)
    {
        const auto n=atlas::make_node(s,i,j);auto r=atlas::trace_cpu(s,n);
        require(n.trace.incident_frame.valid(n.trace.incident_direction),"Bad frame");
        require(r.audit.accepted && r.histogram_relative_error<1e-12,"Node audit");
        auto repeat=atlas::trace_cpu(s,n);require(r.mass_mm2==repeat.mass_mm2,"CPU repeat not bitwise equal");
    }
    s.interior_index=1;
    auto n=atlas::make_node(s,1,1);auto r=atlas::trace_cpu(s,n);
    for(unsigned i=0;i<s.theta_count;++i)for(unsigned j=0;j<s.phi_count;++j)
        require(std::fabs(r.mass_mm2[std::size_t(i)*s.phi_count+j]-(i==0?r.projected_area_mm2/s.phi_count:0))<1e-15,"Index matched delta/cap");
}
void symmetry_test()
{
    auto s=settings();auto n=atlas::make_node(s,1,1);
    for(unsigned rot=0;rot<6;++rot)for(unsigned m=0;m<2;++m)for(unsigned y=0;y<2;++y)
    {
        auto c=n.trace;normalize(transform(c.incident_direction,rot,m,y),c.incident_direction);
        normalize(transform(n.trace.incident_frame.e0,rot,m,y),c.incident_frame.e0);
        c.incident_frame.e1=c.incident_direction.cross(c.incident_frame.e0);
        for(unsigned id: {0u,7u,19u,47u})
        {
            VectorSink a,b;auto e=n.area_rays[id].entry;
            auto aa=trace_from_entry(n.trace,id,e,a);
            e.position=transform(e.position,rot,m,y);
            const auto expected_normal=transform(s.prism.normal(e.face),rot,m,y);
            for(unsigned f=0;f<8;++f)if((expected_normal-s.prism.normal(f)).max_abs()<2e-14){e.face=f;break;}
            auto bb=trace_from_entry(c,id,e,b);
            require(accepted_status(aa.status)&&accepted_status(bb.status)&&a.samples.size()==b.samples.size(),"Symmetry trace status/order");
            for(std::size_t j=0;j<a.samples.size();++j)
            {
                require((transform(a.samples[j].direction,rot,m,y)-b.samples[j].direction).max_abs()<5e-12,"O3 transformed direction");
                require(std::fabs(a.samples[j].power_fraction-b.samples[j].power_fraction)<5e-12,"O3 unpolarized power");
            }
        }
    }
}
void failures_test()
{
    auto s=settings();s.max_internal_hits=1;auto n=atlas::make_node(s,1,1);
    rejects([&]{static_cast<void>(atlas::trace_cpu(s,n));});
    s=settings();s.area_subdivisions=0;rejects([&]{atlas::validate(s);});
    s=settings();s.max_work_bytes=1;rejects([&]{atlas::validate(s);});
}
}
int main()
{
    try{area_test();std::cout<<"area_quadrature PASS\n";nodes_test();std::cout<<"nodes_repeat_delta PASS\n";
        symmetry_test();std::cout<<"all_24_O3_ray_transforms PASS\n";failures_test();std::cout<<"failure_gates PASS\n";return 0;}
    catch(const std::exception& e){std::cerr<<"Atlas host FAILED: "<<e.what()<<'\n';return 1;}
}
