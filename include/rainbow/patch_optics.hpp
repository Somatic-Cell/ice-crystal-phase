#pragma once

#include <rainbow/cuda_module.hpp>
#include <rainbow/device_buffer.hpp>
#include <rainbow/patch_optics_data.hpp>
#include <rainbow/query_direction_grid.hpp>
#include <cstdint>
#include <filesystem>
#include <span>
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
    void evaluate(const PatchAccel& source, const PatchQuery& query,
                  const RaindropTraceConfig& config);
    void write_csv(const std::filesystem::path& path) const;
    [[nodiscard]] bool has_result() const noexcept { return has_result_; }
    [[nodiscard]] const DeviceBuffer<PatchOpticalResult>& results() const noexcept { return results_; }
    [[nodiscard]] std::span<const PatchOpticalResult> host_results() const noexcept { return host_results_; }
    [[nodiscard]] const PatchOpticsStatistics& statistics() const noexcept { return statistics_; }
    void close();
    [[nodiscard]] bool close_noexcept() noexcept;

private:
    void synchronize();
    const CudaContext& cuda_context_;
    CudaModule module_;
    CUfunction function_ = nullptr;
    DeviceBuffer<PatchOpticalResult> results_;
    std::vector<PatchOpticalResult> host_results_;
    std::vector<Vec3> host_directions_; // copy of existing HOST metadata, not GPU readback
    RaindropTraceConfig config_{};
    QueryDirectionGrid grid_{};
    PatchOpticsStatistics statistics_{};
    bool has_grid_ = false, has_result_ = false, has_pending_work_ = false, is_closed_ = false;
};

} // namespace rainbow
