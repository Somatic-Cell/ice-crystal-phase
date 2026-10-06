#pragma once

#include <rainbow/cuda_module.hpp>
#include <rainbow/device_buffer.hpp>
#include <rainbow/patch_optics_data.hpp>
#include <rainbow/wave_optics_data.hpp>
#include <rainbow/folded_patch_data.hpp>
#include <rainbow/query_direction_grid.hpp>
#include <cstdint>
#include <filesystem>
#include <span>
#include <iosfwd>
#include <vector>

namespace rainbow
{
class PatchAccel;
class PatchQuery;

struct PatchOpticsStatistics
{
    std::uint32_t directions = 0;
    std::uint32_t known_hits_complete_directions = 0;
    std::uint32_t pending_directions = 0;
    std::uint32_t error_directions = 0;
    std::uint64_t evaluated_hits = 0, rejected_hits = 0;
    std::uint64_t missing_source_cells = 0, no_outgoing_source_cells = 0;
};

// Owns only the ordinary CUDA module and optical results. Input buffers are
// borrowed during synchronous evaluate(); they are never downloaded/re-uploaded.
// Caller contract: context, vertices, patches and query must correspond to one
// unchanged trace/build/query result and the same config. Rebuild/retrace/query
// invalidates that combination. Address checks do not detect in-place mutation.
class PatchOptics final
{
public:
    explicit PatchOptics(const CudaContext& cuda_context) noexcept;
    ~PatchOptics() noexcept;
    PatchOptics(const PatchOptics&) = delete;
    PatchOptics& operator=(const PatchOptics&) = delete;
    PatchOptics(PatchOptics&&) = delete;
    PatchOptics& operator=(PatchOptics&&) = delete;

    void load_module(const std::filesystem::path& fatbin_path);
    // Applies to the NEXT evaluate call; legacy comparison can disable the
    // extension without modifying geometry or the original query buffers.
    void enable_folded_patches(bool enabled) noexcept { folded_enabled_ = enabled; }
    void set_folded_patch_config(const FoldedPatchConfig& config);
    [[nodiscard]] const FoldedPatchStatistics& folded_statistics() const noexcept { return folded_statistics_; }
    void evaluate(const PatchAccel& source, const PatchQuery& query,
                  const RaindropTraceConfig& config);
    // Explicit new mode; the original evaluate()/write_csv() semantics remain.
    void evaluate_wave(const PatchAccel& source, const PatchQuery& query,
                       const RaindropTraceConfig& config, const WaveOpticsSettings& settings = {});
    void write_wave_csv(const std::filesystem::path& path) const;
    [[nodiscard]] bool has_wave_result() const noexcept { return has_wave_result_; }
    [[nodiscard]] const WaveOpticsStatistics& wave_statistics() const noexcept { return wave_statistics_; }
    [[nodiscard]] std::span<const FocalOpticalResult> host_focal_results() const noexcept { return host_focal_; }
    [[nodiscard]] std::span<const DiffractionResult> host_diffraction_results() const noexcept { return host_diffraction_; }
    void write_csv(const std::filesystem::path& path) const;
    [[nodiscard]] bool has_result() const noexcept { return has_result_; }
    [[nodiscard]] const DeviceBuffer<PatchOpticalResult>& results() const noexcept { return results_; }
    [[nodiscard]] std::span<const PatchOpticalResult> host_results() const noexcept { return host_results_; }
    [[nodiscard]] const PatchOpticsStatistics& statistics() const noexcept { return statistics_; }
    void close();
    [[nodiscard]] bool close_noexcept() noexcept;

private:
    void synchronize();
    void evaluate_impl(const PatchAccel&, const PatchQuery&, const RaindropTraceConfig&, const WaveOpticsSettings*);
    void run_diffraction();
    FoldedPatchView prepare_folded(const PatchAccel&, const FocalPhaseConfig&, bool with_focal);
    void collect_folded_statistics(bool with_focal);
    void write_folded_metadata(std::ostream&) const;
    const CudaContext& cuda_context_;
    CudaModule module_;
    CUfunction function_ = nullptr;
    DeviceBuffer<PatchOpticalResult> results_;
    DeviceBuffer<FocalOpticalResult> focal_results_;
    DeviceBuffer<RainbowTransition> transitions_;
    DeviceBuffer<DiffractionResult> diffraction_results_;
    DeviceBuffer<std::uint32_t> folded_indices_, folded_written_;
    DeviceBuffer<FoldedPatchRecord> folded_records_;
    FoldedPatchConfig folded_config_{}, result_folded_config_{};
    FoldedPatchStatistics folded_statistics_{};
    bool folded_enabled_ = true, result_used_folded_ = false;
    std::vector<FocalOpticalResult> host_focal_;
    std::vector<RainbowTransition> host_transitions_;
    std::vector<DiffractionResult> host_diffraction_;
    WaveOpticsSettings wave_settings_{};
    DiffractionConfig diffraction_config_{};
    WaveOpticsStatistics wave_statistics_{};
    bool has_wave_result_ = false;
    std::vector<PatchOpticalResult> host_results_;
    std::vector<Vec3> host_directions_; // copy of existing HOST metadata, not GPU readback
    RaindropTraceConfig config_{};
    QueryDirectionGrid grid_{};
    PatchOpticsStatistics statistics_{};
    bool has_grid_ = false, has_result_ = false, has_pending_work_ = false, is_closed_ = false;
};

} // namespace rainbow
