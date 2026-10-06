#pragma once

#include <rainbow/cuda_module.hpp>
#include <rainbow/device_buffer.hpp>
#include <rainbow/optix_context.hpp>
#include <rainbow/outgoing_patch.hpp>

#include <array>
#include <cstdint>
#include <filesystem>
#include <vector>

namespace rainbow
{
class RaindropTracer;

struct PatchBuildStatistics
{
    std::array<std::uint64_t,patch_status_count> cells{};
    std::uint32_t logical_cell_count=0;
    std::uint32_t patch_count=0;
    std::size_t gas_byte_size=0;

    [[nodiscard]] std::uint64_t count(const PatchCellStatus status) const noexcept
    { return cells[static_cast<unsigned int>(status)]; }
};

// Owns the constructed patch records, persistent cell-status map, AABBs and GAS.
// Does NOT own vertices, CudaContext, or OptixContext. All three must outlive
// this object; vertices and the optional second-column buffer must remain
// unchanged and alive until all downstream consumers finish.
// The vertex buffer and both contexts must belong to the same CUDA context.
//
// This is geometric construction, NOT an optical query. Incomplete and folded
// cells are recorded explicitly; a successful build is not a complete phase LUT.
class PatchAccel final
{
public:
    PatchAccel(const CudaContext& cuda_context,const OptixContext& optix_context) noexcept;
    ~PatchAccel() noexcept;
    PatchAccel(const PatchAccel&)=delete;
    PatchAccel& operator=(const PatchAccel&)=delete;
    PatchAccel(PatchAccel&&)=delete;
    PatchAccel& operator=(PatchAccel&&)=delete;

    void load_module(const std::filesystem::path& fatbin_path);

    // Blocking boundary: prior work on cuda_context.stream() is completed;
    // classification, stable compaction and GAS build finish before return.
    // Users of a previous GAS on OTHER streams must already be finished.
    void build(const DeviceBuffer<OutgoingVertex>& vertices,const RaindropTraceConfig& config,
               const DeviceBuffer<Field32>* second_input_fields = nullptr);
    // Preferred production entry: carries both response columns without asking
    // the caller to separately label or forward the second-column view.
    void build(const RaindropTracer& tracer);
    [[nodiscard]] bool is_unpolarized() const noexcept { return source_second_fields_ != 0; }
    [[nodiscard]] CUdeviceptr source_second_fields_address() const noexcept { return source_second_fields_; }

    [[nodiscard]] OptixTraversableHandle handle() const noexcept {return traversable_;}
    [[nodiscard]] bool has_result() const noexcept {return has_result_;}
    [[nodiscard]] const DeviceBuffer<OutgoingPatch>& patches() const noexcept {return patches_;}
    [[nodiscard]] const DeviceBuffer<PatchAabb>& aabbs() const noexcept {return aabbs_;}
    [[nodiscard]] const DeviceBuffer<PatchCellStatus>& cell_statuses() const noexcept {return statuses_;}
    [[nodiscard]] const PatchBuildLayout& layout() const noexcept {return layout_;}
    [[nodiscard]] const PatchBuildStatistics& statistics() const noexcept {return statistics_;}
    [[nodiscard]] CUdeviceptr source_vertices_address() const noexcept {return source_vertices_;}

    // Optional inspection I/O. Normal construction does not read back vertices,
    // patch records, AABBs, or the full cell map.
    void write_csv(const std::filesystem::path& path) const;
    void close();
    [[nodiscard]] bool close_noexcept() noexcept;

private:
    void synchronize();
    void build_gas();
    void clear_result();

    const CudaContext& cuda_context_;
    const OptixContext& optix_context_;
    CudaModule module_;
    DeviceBuffer<PatchCellStatus> statuses_;
    DeviceBuffer<PatchBlockSummary> block_summaries_;
    DeviceBuffer<std::uint32_t> block_offsets_;
    DeviceBuffer<OutgoingPatch> patches_;
    DeviceBuffer<PatchAabb> aabbs_;
    DeviceBuffer<std::byte> gas_scratch_;
    DeviceBuffer<std::byte> gas_output_;

    // Async upload source persists through the destructor's completion wait.
    std::vector<PatchBlockSummary> host_summaries_;
    std::vector<std::uint32_t> host_offsets_;
    PatchBuildLayout layout_{};
    PatchBuildStatistics statistics_{};
    OptixTraversableHandle traversable_=0;
    CUdeviceptr source_vertices_=0,source_second_fields_=0;
    bool has_pending_work_=false;
    bool has_result_=false;
    bool is_closed_=false;
};

} // namespace rainbow
