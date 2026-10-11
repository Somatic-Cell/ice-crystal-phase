#include <ice_crystal/atlas.hpp>
#include <rainbow/phase_cdf_math.hpp>
#include <algorithm>
#include <array>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace iceCrystal::atlas
{
std::vector<AreaRay> make_area_rule(const HexPrism& h,Vec3 k,std::uint32_t n)
{
    if(!h.valid()||!is_unit(k)||!n||std::uint64_t(n)*n>1000000)
        throw std::invalid_argument("Invalid area quadrature geometry / subdivision count.");
    const double ga=1/std::sqrt(3.0),gb=std::sqrt(3.0/5.0);
    const std::array<double,2> u{(1-ga)/2,(1+ga)/2};
    const std::array<double,3> v{(1-gb)/2,0.5,(1+gb)/2},vw{5.0/18,4.0/9,5.0/18};
    std::vector<AreaRay> rays;
    Sum area{};
    for(unsigned f=0;f<8;++f)
    {
        const double projected=h.face_area_mm2(f)*max_value(0,-h.normal(f).dot(k));
        if(projected==0)continue;
        const unsigned sectors=f<6?1:6;
        for(unsigned sector=0;sector<sectors;++sector)
        for(std::uint32_t a=0;a<n;++a)for(std::uint32_t b=0;b<n;++b)
        for(unsigned i=0;i<2;++i)for(unsigned j=0;j<3;++j)
        {
            const double x=(a+u[i])/n,y=(b+v[j])/n;
            AreaRay q{};q.entry.face=f;q.entry.valid=true;
            if(f<6)
            {
                const auto normal=h.normal(f);const Vec3 t{-normal.z,0,normal.x};
                q.entry.position=normal*(sqrt3/2)+t*(x-0.5)+Vec3{0,(2*y-1)*h.half_length(),0};
            }
            else
            {
                q.entry.position=(h.vertex(sector)*(1-y)+h.vertex(sector+1)*y)*std::sqrt(x);
                q.entry.position.y=(f==6?1:-1)*h.half_length();
            }
            q.projected_weight_mm2=(projected/sectors/double(n)/n)*0.5*vw[j];
            if(!(q.projected_weight_mm2>0)||!finite_value(q.projected_weight_mm2)||
               !strict_face_point(h,q.entry.position,f))
                throw std::runtime_error("Area quadrature is not representable / lies on ambiguous boundary.");
            rays.push_back(q);area.add(q.projected_weight_mm2);
        }
    }
    if(rays.empty()||std::fabs(area.value/h.projected_area_mm2(k)-1)>1e-12)
        throw std::runtime_error("Projected-area quadrature failed conservation audit.");
    return rays;
}
void validate(const Settings& s)
{
    if(!s.prism.valid()||!(s.interior_index>0&&s.exterior_index>0&&s.wavelength_nm>0)||
       !finite_value(s.interior_index)||!finite_value(s.exterior_index)||!finite_value(s.wavelength_nm)||
       !s.incident_theta_intervals||!s.incident_phi_intervals||
       s.incident_theta_intervals>4096||s.incident_phi_intervals>4096||
       !s.theta_count||!s.phi_count||std::uint64_t(s.theta_count)*s.phi_count>100000000||
       !s.area_subdivisions||std::uint64_t(s.area_subdivisions)*s.area_subdivisions>1000000||
       !s.max_internal_hits||!s.batch_rays||!s.max_work_bytes)
        throw std::invalid_argument("Invalid atlas settings.");
    for(double v:{s.tail_tolerance,s.max_unresolved_fraction,s.max_balance_error})
        if(!finite_value(v)||v<0||v>=1)throw std::invalid_argument("Invalid audit tolerance.");
    if(s.tail_tolerance>s.max_unresolved_fraction)throw std::invalid_argument("Ray tail exceeds global unresolved budget.");
    const auto cells=std::uint64_t(s.theta_count)*s.phi_count;
    if(s.batch_rays>s.max_work_bytes/sizeof(Sum)/cells)
        throw std::length_error("Requested deterministic per-ray grid exceeds --max-work-mib; lower --batch-rays explicitly.");
}
Node make_node(const Settings& s,std::uint32_t i,std::uint32_t j)
{
    validate(s);
    if(i>s.incident_theta_intervals||j>s.incident_phi_intervals)throw std::out_of_range("Incident node index.");
    Node n;n.theta=(pi/2)*i/s.incident_theta_intervals;n.phi=(pi/6)*j/s.incident_phi_intervals;
    const double st=i==0?0:(i==s.incident_theta_intervals?1:std::sin(n.theta));
    const double ct=i==0?1:(i==s.incident_theta_intervals?0:std::cos(n.theta));
    const double sp=j==0?0:(j==s.incident_phi_intervals?0.5:std::sin(n.phi));
    const double cp=j==0?1:(j==s.incident_phi_intervals?sqrt3/2:std::cos(n.phi));
    Vec3 k{st*cp,ct,st*sp};normalize(k,k);
    n.area_rays=make_area_rule(s.prism,k,s.area_subdivisions);
    n.trace=make_trace_settings(s.prism,k,s.exterior_index,s.interior_index,n.area_rays.size(),0,
        s.max_internal_hits,s.tail_tolerance);
    k=n.trace.incident_direction;
    Vec3 e{ct*cp,-st,ct*sp};normalize(e,e);
    n.trace.incident_frame={e,k.cross(e)};
    if(!n.trace.incident_frame.valid(k))throw std::runtime_error("Invalid incident spherical frame.");
    n.grid={nullptr,k,n.trace.incident_frame,s.theta_count,s.phi_count};
    n.u_edges=rainbow::make_phase_u_edges(s.theta_count);return n;
}
void add_audits(Result& r,std::span<const RayAudit> audits,std::span<const AreaRay> rays,double tolerance)
{
    if(audits.size()!=rays.size())throw std::logic_error("Area/audit size mismatch.");
    const auto summary=summarize(audits,tolerance);merge_summary(r.audit,summary);
    for(std::size_t k=0;k<rays.size();++k)
    {
        const auto& a=audits[k];
        if(!accepted_status(a.status)||!finite_value(a.max_interface_balance_error)||
           a.max_interface_balance_error>tolerance)
        {
            std::ostringstream os;os<<"Deterministic area ray rejected: local_id="<<k<<" face="<<rays[k].entry.face
                <<" status="<<status_name(a.status)<<" residual="<<a.unresolved_power
                <<". No jitter, dropped ray, or renormalized omission is applied.";
            throw std::runtime_error(os.str());
        }
        r.escaped_accumulator.add(rays[k].projected_weight_mm2*a.escaped_power);
        r.unresolved_accumulator.add(rays[k].projected_weight_mm2*a.unresolved_power);
        r.escaped_mm2=r.escaped_accumulator.value;
        r.unresolved_mm2=r.unresolved_accumulator.value;
    }
    require_accepted(summary);
}
void finish(Result& r,const Settings& s,const Node& n)
{
    r.projected_area_mm2=s.prism.projected_area_mm2(n.trace.incident_direction);
    Sum total{};for(double w:r.mass_mm2)
    {if(!finite_value(w)||w<0)throw std::runtime_error("Invalid atlas cell mass.");total.add(w);}
    if(!(total.value>0)||!(r.escaped_mm2>0)||r.audit.rays!=n.area_rays.size())
        throw std::runtime_error("Incomplete atlas node.");
    r.histogram_relative_error=std::fabs(total.value/r.escaped_mm2-1);
    if(r.histogram_relative_error>s.max_balance_error||
       std::fabs((r.escaped_mm2+r.unresolved_mm2)/r.projected_area_mm2-1)>s.max_balance_error||
       r.unresolved_mm2/r.projected_area_mm2>s.max_unresolved_fraction)
        throw std::runtime_error("Atlas area / escaped / residual conservation audit failed.");
}
void validate_node(const Settings& s,const Node& n)
{
    validate(s);
    if(n.grid.nt!=s.theta_count||n.grid.np!=s.phi_count||n.u_edges.size()!=std::size_t(s.theta_count)+1||
       n.area_rays.empty()||n.trace.total_incident_samples!=n.area_rays.size()||
       n.trace.max_internal_hits!=s.max_internal_hits||n.trace.residual_power_tolerance!=s.tail_tolerance||
       n.trace.prism.circumradius_mm!=s.prism.circumradius_mm||n.trace.prism.length_mm!=s.prism.length_mm||
       n.trace.interior_index!=s.interior_index||n.trace.exterior_index!=s.exterior_index||
       !n.grid.frame.valid(n.grid.ki)||!settings_valid(n.trace)||
       (n.grid.ki-n.trace.incident_direction).max_abs()!=0)
        throw std::invalid_argument("Atlas node does not match trace settings.");
    if(n.u_edges.front()!=0||n.u_edges.back()!=1)throw std::invalid_argument("Invalid u endpoints.");
    for(std::size_t i=1;i<n.u_edges.size();++i)
        if(!finite_value(n.u_edges[i])||!(n.u_edges[i]>n.u_edges[i-1]))throw std::invalid_argument("Invalid u edges.");
}
Result trace_cpu(const Settings& s,const Node& n)
{
    validate_node(s,n);
    Result r;const auto size=std::size_t(s.theta_count)*s.phi_count;
    std::vector<Sum> accumulated(size),per_ray(size);
    for(std::size_t id=0;id<n.area_rays.size();++id)
    {
        std::fill(per_ray.begin(),per_ray.end(),Sum{});
        MassSink sink{n.view(),per_ray.data(),1,0,n.area_rays[id].projected_weight_mm2};
        const auto audit=trace_from_entry(n.trace,id,n.area_rays[id].entry,sink);
        add_audits(r,{&audit,1},{&n.area_rays[id],1},s.max_balance_error);
        for(std::size_t k=0;k<size;++k)if(per_ray[k].value>0)accumulated[k].add(per_ray[k].value);
    }
    r.mass_mm2.resize(size);for(std::size_t k=0;k<size;++k)r.mass_mm2[k]=accumulated[k].value;
    finish(r,s,n);return r;
}
}
