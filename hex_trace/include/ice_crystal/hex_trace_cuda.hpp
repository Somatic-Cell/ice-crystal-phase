#pragma once
#include <ice_crystal/hex_trace.hpp>
#include <rainbow/cuda_driver.hpp>
#include <rainbow/cuda_module.hpp>
#include <rainbow/device_buffer.hpp>
#include <filesystem>

namespace iceCrystal
{
// No CUDA Runtime API. Owns output storage; context is borrowed and must outlive
// this object. Each trace() is synchronous. No automatic CPU fallback.
class HexPrismTracer final
{
public:
    explicit HexPrismTracer(const rainbow::CudaContext&);
    ~HexPrismTracer() noexcept;
    HexPrismTracer(const HexPrismTracer&)=delete;
    HexPrismTracer& operator=(const HexPrismTracer&)=delete;
    void load_module(const std::filesystem::path& fatbin);
    void trace(const TraceSettings&,std::uint64_t first,std::size_t count,
               std::size_t maximum_output_records=4*1024*1024);
    [[nodiscard]] TraceBatch download() const;
    [[nodiscard]] std::vector<RayAudit> download_audits() const;
    // Device-resident CSR output for a later phase accumulator. Must not outlive
    // the tracer or the next call to trace(). Addresses alone do not certify quality.
    [[nodiscard]] const rainbow::DeviceBuffer<OutgoingSample>& outgoing() const noexcept {return outgoing_;}
    [[nodiscard]] const rainbow::DeviceBuffer<std::uint64_t>& offsets() const noexcept {return offsets_;}
    [[nodiscard]] const rainbow::DeviceBuffer<RayAudit>& audits() const noexcept {return audits_;}
    [[nodiscard]] bool has_result() const noexcept {return has_result_;}
private:
    void synchronize();
    void launch(const char*,void**,unsigned blocks);
    const rainbow::CudaContext& context_;
    rainbow::CudaModule module_;
    rainbow::DeviceBuffer<RayAudit> audits_;
    rainbow::DeviceBuffer<std::uint64_t> offsets_,total_;
    rainbow::DeviceBuffer<OutgoingSample> outgoing_;
    std::size_t count_=0;
    std::uint64_t first_=0;
    bool pending_=false,has_result_=false,module_ready_=false;
};
} // namespace iceCrystal
