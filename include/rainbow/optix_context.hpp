#pragma once

#include <rainbow/cuda_driver.hpp>
#include <optix.h>

namespace rainbow
{
// CUDA context は借用，OptixDeviceContext は所有する．stream / primary retain は解放しない．
class OptixContext final
{
public:
    explicit OptixContext(const CudaContext& cuda_context);
    ~OptixContext() noexcept;
    OptixContext(const OptixContext&)=delete;
    OptixContext& operator=(const OptixContext&)=delete;
    OptixContext(OptixContext&&)=delete;
    OptixContext& operator=(OptixContext&&)=delete;

    // 旧 smoke と本番追跡で，function table の初期化を共有する．
    static void initialize_runtime();
    [[nodiscard]] OptixDeviceContext handle() const noexcept {return context_;}
    [[nodiscard]] bool close_noexcept() noexcept;
private:
    static void log_callback(unsigned level,const char* tag,const char* message,void*) noexcept;
    const CudaContext& cuda_context_;
    OptixDeviceContext context_=nullptr;
};
}
