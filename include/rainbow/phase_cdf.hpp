#pragma once

#include <rainbow/phase_cdf_data.hpp>
#include <rainbow/cuda_module.hpp>
#include <rainbow/device_buffer.hpp>
#include <filesystem>
#include <vector>

namespace rainbow
{
class PatchOptics;
class RaindropTracer;

// Driver-API-only owner. Optical input buffers are borrowed during synchronous
// build(); no raw optical result is downloaded or passed through a CSV file.
class PhaseCdf final
{
public:
    explicit PhaseCdf(const CudaContext& context);
    ~PhaseCdf() noexcept;
    PhaseCdf(const PhaseCdf&) = delete;
    PhaseCdf& operator=(const PhaseCdf&) = delete;
    void load_module(const std::filesystem::path& path);
    void validate_trace(const RaindropTracer& trace);
    void build(const PatchOptics&, std::uint32_t ntheta, std::uint32_t nphi,
               PhaseDensityStage stage = PhaseDensityStage::Diffraction,
               const PhaseCdfPolicy& policy = {});
    void build(PhaseDensityView input, std::uint32_t ntheta, std::uint32_t nphi,
               const PhaseCdfPolicy& policy = {});
    // Writes only the three standard NPY files; metadata/commit are the caller's responsibility.
    void write_arrays(const std::filesystem::path& staging_directory) const;
    [[nodiscard]] const PhaseCdfReduction& input_statistics() const noexcept { return input_statistics_; }
    [[nodiscard]] const PhaseCdfReduction& audit_statistics() const noexcept { return audit_statistics_; }
    [[nodiscard]] const std::vector<double>& u_edges() const noexcept { return host_edges_; }
    [[nodiscard]] const DeviceBuffer<double>& phi_cdf() const noexcept { return phi_; }
    [[nodiscard]] const DeviceBuffer<double>& theta_cdf() const noexcept { return theta_; }
    [[nodiscard]] bool valid() const noexcept { return valid_; }
private:
    void synchronize();
    void allocate_partials(std::uint64_t count);
    void reduce_to(DeviceBuffer<PhaseCdfReduction>& output);
    void launch(const char* name, void** args, unsigned gx, unsigned gy=1,
                unsigned tx=256, unsigned ty=1);
    const CudaContext& context_;
    CudaModule module_;
    DeviceBuffer<double> edges_, phi_, theta_, column_sums_, total_;
    DeviceBuffer<PhaseCdfReduction> partials_, input_summary_, audit_summary_;
    std::vector<double> host_edges_;
    PhaseCdfReduction input_statistics_{}, audit_statistics_{};
    std::uint32_t nt_=0, np_=0, partial_count_=0;
    bool pending_=false, valid_=false;
};
} // namespace rainbow
