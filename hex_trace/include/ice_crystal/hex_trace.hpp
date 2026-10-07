#pragma once
#include <ice_crystal/hex_trace_core.hpp>
#include <cstddef>
#include <span>
#include <vector>

namespace iceCrystal
{
// Host validation and CPU reference use the SAME core as the CUDA kernels.
// Direction normalization is explicit here; no hidden fallback direction.
TraceSettings make_trace_settings(HexPrism prism,Vec3 incident,double exterior_index,
    double interior_index,std::uint64_t total_samples,std::uint64_t seed,
    std::uint32_t max_internal_hits,double residual_power_tolerance);
void validate_settings(const TraceSettings&);
void validate_batch_range(const TraceSettings&,std::uint64_t first,std::size_t count);
const char* status_name(TraceStatus) noexcept;
struct TraceBatch
{
    std::uint64_t first_sample=0;
    std::vector<RayAudit> audits;
    std::vector<std::uint64_t> offsets; // CSR range [offsets[i], offsets[i+1])
    std::vector<OutgoingSample> outgoing;
};
TraceBatch trace_batch_cpu(const TraceSettings&,std::uint64_t first,std::size_t count,
    std::size_t maximum_output_records=4*1024*1024);
struct TraceSummary
{
    std::uint64_t rays=0,outputs=0,status_counts[8]{};
    double escaped_power_sum=0,unresolved_power_sum=0;
    double maximum_balance_error=0,maximum_interface_balance_error=0;
    bool accepted=true;
};
TraceSummary summarize(std::span<const RayAudit> audits,double balance_tolerance=1e-10);
void merge_summary(TraceSummary& destination,const TraceSummary& source);
void require_accepted(const TraceSummary&);
} // namespace iceCrystal
