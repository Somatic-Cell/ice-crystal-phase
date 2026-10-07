#pragma once
#include <ice_crystal/hex_prism.hpp>
#include <ice_crystal/flux_interface.hpp>
#include <type_traits>

namespace iceCrystal
{
enum class TraceStatus : std::uint32_t
{
    Exhausted=0, TailToleranceReached=1, InteractionLimit=2,
    AmbiguousBoundary=3, InvalidInput=4, NumericalFailure=5,
    OutputOverflow=6, ReplayMismatch=7
};
enum class ExitKind : std::uint32_t { ExteriorReflection=0, InteriorTransmission=1 };
struct TraceSettings
{
    HexPrism prism{};
    Vec3 incident_direction{};
    Frame incident_frame{};
    // Evaluated indices supplied by the caller for its wavelength. No built-in
    // constant is silently presented as ice dispersion or birefringent ice.
    double exterior_index=0,interior_index=0;
    double residual_power_tolerance=1e-12; // absolute fraction of each unit ray
    std::uint64_t seed=0,total_incident_samples=0;
    std::uint32_t max_internal_hits=256;
    std::uint32_t reserved=0;
};
struct OutgoingSample
{
    Vec3 position_mm,direction;
    Frame frame;
    JonesFlux response;
    double power_fraction=0;
    double internal_length_mm=0;
    std::uint64_t incident_sample_id=0;
    std::uint32_t entry_face=8,exit_face=8,internal_reflections=0;
    ExitKind kind=ExitKind::ExteriorReflection;
    // Physical cell contribution (before orientation weighting):
    // projected_area_mm2 / total_incident_samples * power_fraction.
};
struct RayAudit
{
    double escaped_power=0,unresolved_power=1,balance_error=0;
    double max_interface_balance_error=0;
    std::uint64_t output_count=0;
    std::uint32_t internal_hits=0,tir_count=0;
    TraceStatus status=TraceStatus::InvalidInput;
    std::uint32_t reserved=0;
};
struct CountSink
{
    ICE_HD bool push(const OutgoingSample&) noexcept {return true;}
};
struct BufferSink
{
    OutgoingSample* data=nullptr;
    std::uint64_t capacity=0,size=0;
    ICE_HD bool push(const OutgoingSample& s) noexcept
    {
        if(size>=capacity||!data)return false;
        data[size++]=s;return true;
    }
};
ICE_HD inline bool accepted_status(TraceStatus s) noexcept
{ return s==TraceStatus::Exhausted || s==TraceStatus::TailToleranceReached; }
ICE_HD inline bool settings_valid(const TraceSettings& c) noexcept
{
    const double a=c.prism.projected_area_mm2(c.incident_direction);
    return c.prism.valid()&&c.incident_frame.valid(c.incident_direction)&&
        finite_value(c.exterior_index)&&finite_value(c.interior_index)&&c.exterior_index>0&&c.interior_index>0&&
        finite_value(c.residual_power_tolerance)&&c.residual_power_tolerance>=0&&c.residual_power_tolerance<1&&
        c.total_incident_samples>0&&c.max_internal_hits>0&&finite_value(a)&&a>0;
}
ICE_HD inline RayAudit finish_audit(RayAudit a,const Sum& escaped,double residual,TraceStatus status) noexcept
{
    a.escaped_power=escaped.value;a.unresolved_power=residual;
    a.balance_error=escaped.value+residual-1;a.status=status;return a;
}

template<class Sink>
ICE_HD TraceStatus emit_exit(Sink& sink,RayAudit& audit,Sum& escaped,const TraceSettings& c,
    std::uint64_t id,const FluxRay& ray,Vec3 position,unsigned entry,unsigned exit,
    unsigned reflections,ExitKind kind,double length) noexcept
{
    const double power=ray.response.power();
    if(!finite_value(power)||power<0||(!power && ray.response.nonzero()))return TraceStatus::NumericalFailure;
    if(power==0)return TraceStatus::Exhausted;
    OutgoingSample s{};
    s.position_mm=position*c.prism.circumradius_mm;s.direction=ray.direction;s.frame=ray.frame;
    s.response=ray.response;s.power_fraction=power;s.internal_length_mm=length;
    s.incident_sample_id=id;s.entry_face=entry;s.exit_face=exit;s.internal_reflections=reflections;s.kind=kind;
    if(!s.position_mm.is_finite())return TraceStatus::NumericalFailure;
    if(!sink.push(s))return TraceStatus::OutputOverflow;
    escaped.add(power);++audit.output_count;return TraceStatus::Exhausted;
}

// Deterministic splitting, NOT a single randomly chosen path. For one convex
// prism, exterior branches cannot re-enter it: emit exterior R once, then emit
// T at every internal encounter and continue only internal R. Thus all retained
// orders require a LINEAR chain, not an exponentially growing tree.
template<class Sink>
ICE_HD RayAudit trace_from_entry(const TraceSettings& c,std::uint64_t id,EntrySample entry,Sink& sink) noexcept
{
    RayAudit a{};Sum escaped{};
    if(!settings_valid(c)||id>=c.total_incident_samples||!entry.valid||entry.face>=8)
        return finish_audit(a,escaped,1,TraceStatus::InvalidInput);
    if(!strict_face_point(c.prism,entry.position,entry.face))
        return finish_audit(a,escaped,1,TraceStatus::AmbiguousBoundary);
    if(!(c.prism.normal(entry.face).dot(c.incident_direction)<0))
        return finish_audit(a,escaped,1,TraceStatus::InvalidInput);
    FluxRay incident{c.incident_direction,c.incident_frame,JonesFlux::identity()};
    auto split=split_flux(incident,c.prism.normal(entry.face),c.exterior_index,c.interior_index);
    if(!split.valid)return finish_audit(a,escaped,1,TraceStatus::NumericalFailure);
    a.max_interface_balance_error=::fabs(split.balance_error);
    const auto entry_emission=emit_exit(sink,a,escaped,c,id,split.reflected,entry.position,entry.face,entry.face,
                  0,ExitKind::ExteriorReflection,0);
    if(entry_emission!=TraceStatus::Exhausted)
        return finish_audit(a,escaped,1,entry_emission);
    if(!split.has_transmission)
    {
        a.tir_count=1;
        return finish_audit(a,escaped,0,TraceStatus::Exhausted);
    }
    FluxRay ray=split.transmitted;
    Vec3 p=entry.position;
    unsigned current_face=entry.face,reflections=0;
    Sum length{};
    for(unsigned encounter=0;encounter<c.max_internal_hits;++encounter)
    {
        const double power=ray.response.power();
        if(!finite_value(power)||power<0||(!power&&ray.response.nonzero()))
            return finish_audit(a,escaped,power,TraceStatus::NumericalFailure);
        if(power==0)return finish_audit(a,escaped,0,TraceStatus::Exhausted);
        if(power<=c.residual_power_tolerance)
            return finish_audit(a,escaped,power,TraceStatus::TailToleranceReached);
        auto hit=next_boundary(c.prism,p,ray.direction,current_face);
        if(hit.status!=HitStatus::Regular)
            return finish_audit(a,escaped,power,hit.status==HitStatus::AmbiguousBoundary?
                TraceStatus::AmbiguousBoundary:TraceStatus::NumericalFailure);
        const double segment=hit.distance*c.prism.circumradius_mm;
        if(!finite_value(segment)||segment<0)
            return finish_audit(a,escaped,power,TraceStatus::NumericalFailure);
        length.add(segment);
        if(!finite_value(length.value))return finish_audit(a,escaped,power,TraceStatus::NumericalFailure);
        split=split_flux(ray,-c.prism.normal(hit.face),c.interior_index,c.exterior_index);
        if(!split.valid)return finish_audit(a,escaped,power,TraceStatus::NumericalFailure);
        ++a.internal_hits;
        if(split.total_internal_reflection)++a.tir_count;
        a.max_interface_balance_error=max_value(a.max_interface_balance_error,::fabs(split.balance_error));
        if(split.has_transmission)
        {
            const auto emission=emit_exit(sink,a,escaped,c,id,split.transmitted,hit.position,
                entry.face,hit.face,reflections,ExitKind::InteriorTransmission,length.value);
            if(emission!=TraceStatus::Exhausted)return finish_audit(a,escaped,power,emission);
        }
        ray=split.reflected;p=hit.position;current_face=hit.face;++reflections;
    }
    const double remaining=ray.response.power();
    if(!finite_value(remaining)||remaining<0||(!remaining&&ray.response.nonzero()))
        return finish_audit(a,escaped,remaining,TraceStatus::NumericalFailure);
    return finish_audit(a,escaped,remaining,remaining==0?TraceStatus::Exhausted:
        (remaining<=c.residual_power_tolerance?TraceStatus::TailToleranceReached:TraceStatus::InteractionLimit));
}
template<class Sink>
ICE_HD RayAudit trace_sample(const TraceSettings& c,std::uint64_t id,Sink& sink) noexcept
{
    const auto entry=sample_entry(c.prism,c.incident_direction,
        sample_uniform(c.seed,id,0),sample_uniform(c.seed,id,1),sample_uniform(c.seed,id,2));
    return trace_from_entry(c,id,entry,sink);
}
ICE_HD inline bool replay_equal(const RayAudit& a,const RayAudit& b) noexcept
{
    return a.output_count==b.output_count && a.status==b.status &&
        a.internal_hits==b.internal_hits && a.tir_count==b.tir_count &&
        ::fabs(a.escaped_power-b.escaped_power)<=1e-13 &&
        ::fabs(a.unresolved_power-b.unresolved_power)<=1e-13;
}
static_assert(std::is_trivially_copyable_v<TraceSettings>);
static_assert(std::is_trivially_copyable_v<OutgoingSample>);
static_assert(std::is_trivially_copyable_v<RayAudit>);
} // namespace iceCrystal
