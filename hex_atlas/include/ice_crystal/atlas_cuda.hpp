#pragma once
#include <ice_crystal/atlas.hpp>
#include <rainbow/cuda_module.hpp>
#include <rainbow/device_buffer.hpp>
namespace iceCrystal::atlas
{
// Fixed per-ray columns + one writer per output bin: no floating atomicAdd.
class CudaTracer final
{
public:
    CudaTracer(const rainbow::CudaContext&,const std::filesystem::path& module);
    ~CudaTracer() noexcept;
    CudaTracer(const CudaTracer&)=delete;
    CudaTracer& operator=(const CudaTracer&)=delete;
    Result trace(const Settings&,const Node&);
private:
    void sync();
    void launch(const char*,void**,unsigned);
    const rainbow::CudaContext& context_;
    rainbow::CudaModule module_;
    rainbow::DeviceBuffer<AreaRay> rays_;
    rainbow::DeviceBuffer<RayAudit> audits_;
    rainbow::DeviceBuffer<Sum> partial_,mass_;
    rainbow::DeviceBuffer<double> edges_;
    bool pending_=false;
};
}
