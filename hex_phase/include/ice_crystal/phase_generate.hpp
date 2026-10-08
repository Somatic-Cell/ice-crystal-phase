#pragma once
#include <ice_crystal/orientation_plan.hpp>
#include <ice_crystal/phase_table.hpp>
#include <ice_crystal/hex_trace.hpp>
#include <functional>

namespace iceCrystal
{
struct PhaseGenerationSettings
{
    HexPrism prism{};
    PhaseGrid grid{};
    double wavelength_nm=0,interior_index=0,exterior_index=1;
    std::uint32_t max_internal_hits=1024,out_nt=0,out_np=0;
    double tail_tolerance=1e-12,max_unresolved_fraction=1e-10;
    double max_balance_error=1e-10,max_histogram_relative_error=1e-10;
    std::size_t batch_size=4096,max_output_records=4*1024*1024;
    CdfPolicy cdf_policy{};
};
struct CompositionAudit
{
    TraceSummary trace{};
    std::uint64_t orientations_done=0;
    Sum expected_area,escaped_area,unresolved_area;
    HistogramStats histogram{};
    HistogramMoments grid_moments{};
    double histogram_relative_error=0;
};
void validate_generation_settings(const PhaseGenerationSettings&);
TraceSettings node_trace_settings(const PhaseGenerationSettings&,const OrientationNode&);
double node_area(const PhaseGenerationSettings&,const OrientationNode&);
TraceSummary summarize_phase_rays(std::span<const RayAudit>,double balance_tolerance=1e-10);
void record_batch_audit(CompositionAudit&,const TraceSummary&,double scale);
void validate_composition(CompositionAudit&,const PhaseGenerationSettings&,const OrientationPlan&);
void write_phase_metadata(const std::filesystem::path&,const PhaseGenerationSettings&,
    const OrientationPlan&,const CompositionAudit&,const CdfReport&,const char* backend);
void write_failure(const std::filesystem::path&,std::string_view,const CompositionAudit&);
} // namespace iceCrystal
