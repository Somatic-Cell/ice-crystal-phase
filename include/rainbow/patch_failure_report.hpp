#pragma once

#include <rainbow/patch_regularity_audit.hpp>
#include <rainbow/wave_optics_data.hpp>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

namespace rainbow
{
class CudaContext;
class PatchAccel;
class PatchQuery;
class PatchOptics;

// Host-owned witness records, not kernel launch structures.
struct FailedDirectionWitness
{
    std::uint32_t direction_id=0;
    Vec3 direction{};
    PatchQuerySummary query{};
    PatchOpticalResult optical{};
    FocalOpticalResult focal{};
    DiffractionResult diffraction{};
    std::vector<PatchQueryHit> hits; // ALL valid hit records of this direction
};
struct FailedPatchWitness
{
    std::uint32_t compact_index=0;
    OutgoingPatch patch{}; // original GPU classification is authoritative
    std::array<OutgoingVertex,4> vertices{};
    std::array<Field32,4> second_input_fields{};
};
enum class PatchWitnessOrigin { SyntheticTest, GpuCapture, CpuRetrace };
struct PatchFailureSnapshot
{
    PatchWitnessOrigin origin = PatchWitnessOrigin::SyntheticTest;
    IncidentPolarization input_polarization = IncidentPolarization::SingleJones;
    RaindropTraceConfig config{};
    WaveOpticsSettings wave_settings{};
    std::uint32_t theta_count=0,phi_count=0,stored_patch_count=0;
    std::uint64_t missing_source_cells=0,no_outgoing_source_cells=0;
    std::uint64_t query_error_directions=0,optical_incomplete_directions=0;
    std::uint64_t focal_unavailable_directions=0,diffraction_unavailable_directions=0;
    std::uint64_t nonfinite_directions=0,underresolved_directions=0;
    std::uint64_t gpu_read_bytes=0,gpu_read_calls=0;
    std::vector<FailedDirectionWitness> directions;
    std::vector<FailedPatchWitness> patches;
    // first_problem IDs without a compact record (e.g. a corrupt query record).
    std::vector<std::uint32_t> absent_problem_patch_ids;
};

// Captures only witnesses used by failed/pending directions. Owns the host copy,
// not CUDA/OptiX resources. Does NOT change classification, optical values or flags.
class PatchFailureReport final
{
public:
    explicit PatchFailureReport(PatchFailureSnapshot snapshot);
    [[nodiscard]] static PatchFailureReport capture(
        const CudaContext& context,const PatchAccel& patches,const PatchQuery& query,
        const PatchOptics& optics,const RaindropTraceConfig& config,
        const WaveOpticsSettings& wave_settings);
    void write_json(const std::filesystem::path& path) const;
    [[nodiscard]] const PatchFailureSnapshot& snapshot() const noexcept {return snapshot_;}
private:
    PatchFailureSnapshot snapshot_;
};
} // namespace rainbow
