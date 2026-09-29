#include <rainbow/optix_context.hpp>
#include <rainbow/cuda_error.hpp>
#include <rainbow/optix_error.hpp>

// 定義はこの translation unit だけ．既存 src/optix.cpp からは削除する．
#include <optix_function_table_definition.h>
#include <optix_stubs.h>
#include <cstdio>
#include <mutex>

namespace rainbow
{
void OptixContext::initialize_runtime()
{
    static std::once_flag flag;
    std::call_once(flag,[]{RAINBOW_OPTIX_CHECK(optixInit());});
}
void OptixContext::log_callback(unsigned level,const char* tag,const char* message,void*) noexcept
{
    std::fprintf(stderr,"[OptiX][%u][%s] %s\n",level,tag?tag:"",message?message:"");
}
OptixContext::OptixContext(const CudaContext& cuda_context):cuda_context_(cuda_context)
{
    cuda_context_.make_current();
    initialize_runtime();
    OptixDeviceContextOptions options{};
    options.logCallbackFunction=&log_callback;
    options.logCallbackLevel=4;
#ifndef NDEBUG
    options.validationMode=OPTIX_DEVICE_CONTEXT_VALIDATION_MODE_ALL;
#endif
    RAINBOW_OPTIX_CHECK(optixDeviceContextCreate(cuda_context_.handle(),&options,&context_));
}
OptixContext::~OptixContext() noexcept {static_cast<void>(close_noexcept());}
bool OptixContext::close_noexcept() noexcept
{
    if(!context_) return true;
    if(!detail::report_cuda_cleanup_result(cuCtxSetCurrent(cuda_context_.handle()),"cuCtxSetCurrent(OptixContext)")) return false;
    const auto result=optixDeviceContextDestroy(context_);
    context_=nullptr;
    return detail::report_optix_cleanup_result(result,"optixDeviceContextDestroy");
}
}
