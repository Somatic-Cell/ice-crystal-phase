#include <ice_crystal/hex_trace.hpp>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>

namespace iceCrystal
{
void validate_settings(const TraceSettings& s)
{
    if(!settings_valid(s))throw std::invalid_argument("Invalid hex trace geometry, direction/frame, indices, sample count or termination policy.");
    const double area=s.prism.projected_area_mm2(s.incident_direction);
    if(!(area/double(s.total_incident_samples)>0))
        throw std::invalid_argument("Per-incident-sample area underflows FP64.");
}
TraceSettings make_trace_settings(HexPrism p,Vec3 incident,double ni,double nt,
    std::uint64_t total,std::uint64_t seed,std::uint32_t max_hits,double tail)
{
    TraceSettings s{};s.prism=p;
    if(!normalize(incident,s.incident_direction)||!make_frame(s.incident_direction,s.incident_frame))
        throw std::invalid_argument("Incident direction must be finite and nonzero.");
    s.exterior_index=ni;s.interior_index=nt;s.total_incident_samples=total;s.seed=seed;
    s.max_internal_hits=max_hits;s.residual_power_tolerance=tail;validate_settings(s);return s;
}
void validate_batch_range(const TraceSettings& s,std::uint64_t first,std::size_t n)
{
    validate_settings(s);
    if(n==0||first>=s.total_incident_samples||n>s.total_incident_samples-first)
        throw std::invalid_argument("Batch is outside the global incident-sample range.");
    if(n==(std::numeric_limits<std::size_t>::max)())throw std::length_error("Batch offset-array size overflow.");
}
const char* status_name(TraceStatus s) noexcept
{
    switch(s)
    {
    case TraceStatus::Exhausted:return "exhausted";
    case TraceStatus::TailToleranceReached:return "tail_tolerance_reached";
    case TraceStatus::InteractionLimit:return "interaction_limit";
    case TraceStatus::AmbiguousBoundary:return "ambiguous_boundary";
    case TraceStatus::InvalidInput:return "invalid_input";
    case TraceStatus::NumericalFailure:return "numerical_failure";
    case TraceStatus::OutputOverflow:return "output_overflow";
    case TraceStatus::ReplayMismatch:return "replay_mismatch";
    }
    return "unknown";
}
TraceBatch trace_batch_cpu(const TraceSettings& s,std::uint64_t first,std::size_t n,std::size_t max_records)
{
    validate_batch_range(s,first,n);
    TraceBatch b{};b.first_sample=first;b.audits.resize(n);b.offsets.resize(n+1);
    for(std::size_t i=0;i<n;++i)
    {
        CountSink sink{};b.audits[i]=trace_sample(s,first+i,sink);
        const auto count=b.audits[i].output_count;
        if(count>max_records-b.offsets[i])throw std::length_error("Output memory budget exceeded; batch not silently reduced.");
        b.offsets[i+1]=b.offsets[i]+count;
    }
    if(b.offsets.back()>(std::numeric_limits<std::size_t>::max)()/sizeof(OutgoingSample))
        throw std::length_error("Outgoing byte count overflow.");
    b.outgoing.resize(static_cast<std::size_t>(b.offsets.back()));
    for(std::size_t i=0;i<n;++i)
    {
        auto count=b.offsets[i+1]-b.offsets[i];
        auto* data=count?b.outgoing.data()+b.offsets[i]:nullptr;
        BufferSink sink{data,count,0};
        const auto audit=trace_sample(s,first+i,sink);
        if(!replay_equal(b.audits[i],audit)||sink.size!=count)
            b.audits[i].status=TraceStatus::ReplayMismatch;
    }
    return b;
}
TraceSummary summarize(std::span<const RayAudit> audits,double tolerance)
{
    if(!finite_value(tolerance)||tolerance<0)throw std::invalid_argument("Invalid balance tolerance.");
    TraceSummary s{};Sum escaped{},tail{};
    for(const auto& a:audits)
    {
        const auto status=static_cast<unsigned>(a.status);
        if(status>=8)throw std::runtime_error("Unknown audit status / incompatible device module.");
        ++s.rays;s.outputs+=a.output_count;++s.status_counts[status];
        escaped.add(a.escaped_power);tail.add(a.unresolved_power);
        s.maximum_balance_error=(std::max)(s.maximum_balance_error,::fabs(a.balance_error));
        s.maximum_interface_balance_error=(std::max)(s.maximum_interface_balance_error,a.max_interface_balance_error);
        s.accepted=s.accepted && accepted_status(a.status) && finite_value(a.escaped_power)&&finite_value(a.unresolved_power)&&
            a.escaped_power>=0&&a.unresolved_power>=0&&finite_value(a.balance_error)&&::fabs(a.balance_error)<=tolerance;
    }
    s.escaped_power_sum=escaped.value;s.unresolved_power_sum=tail.value;return s;
}
void merge_summary(TraceSummary& d,const TraceSummary& s)
{
    if(s.rays>(std::numeric_limits<std::uint64_t>::max)()-d.rays ||
       s.outputs>(std::numeric_limits<std::uint64_t>::max)()-d.outputs)
        throw std::overflow_error("Summary count overflow.");
    d.rays+=s.rays;d.outputs+=s.outputs;
    for(unsigned i=0;i<8;++i)d.status_counts[i]+=s.status_counts[i];
    d.escaped_power_sum+=s.escaped_power_sum;d.unresolved_power_sum+=s.unresolved_power_sum;
    d.maximum_balance_error=(std::max)(d.maximum_balance_error,s.maximum_balance_error);
    d.maximum_interface_balance_error=(std::max)(d.maximum_interface_balance_error,s.maximum_interface_balance_error);
    d.accepted=d.accepted&&s.accepted;
}
void require_accepted(const TraceSummary& s)
{
    if(!s.rays||!s.accepted)throw std::runtime_error("Hex trace is NOT accepted. Inspect per-ray statuses and unresolved power; do not renormalize it away.");
}
} // namespace iceCrystal
