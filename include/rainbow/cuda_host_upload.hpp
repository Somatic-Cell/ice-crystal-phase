#pragma once

#include <rainbow/cuda_driver.hpp>
#include <rainbow/cuda_error.hpp>
#include <cstddef>
#include <cstring>
#include <stdexcept>

namespace rainbow
{
namespace detail
{
// Owns page-locked staging until the explicit consumer stream has completed.
class PinnedHostUpload final
{
public:
    PinnedHostUpload(const CudaContext& context, std::size_t bytes)
        : context_(context)
    {
        context_.make_current();
        RAINBOW_CUDA_CHECK(cuMemHostAlloc(&data_, bytes, 0));
    }
    ~PinnedHostUpload() noexcept
    {
        if(pending_)
        {
            static_cast<void>(report_cuda_cleanup_result(
                cuCtxSetCurrent(context_.handle()), "upload cleanup context"));
            static_cast<void>(report_cuda_cleanup_result(
                cuStreamSynchronize(context_.stream()), "upload cleanup stream"));
        }
        if(data_) static_cast<void>(report_cuda_cleanup_result(
            cuMemFreeHost(data_), "cuMemFreeHost(upload staging)"));
    }
    PinnedHostUpload(const PinnedHostUpload&) = delete;
    PinnedHostUpload& operator=(const PinnedHostUpload&) = delete;

    void copy_to(CUdeviceptr destination, const void* source, std::size_t bytes)
    {
        std::memcpy(data_, source, bytes);
        pending_ = true;
        RAINBOW_CUDA_CHECK(cuMemcpyHtoDAsync(
            destination, data_, bytes, context_.stream()));
        RAINBOW_CUDA_CHECK(cuStreamSynchronize(context_.stream()));
        pending_ = false;
    }
private:
    const CudaContext& context_;
    void* data_ = nullptr;
    bool pending_ = false;
};
} // namespace detail

// Synchronous with respect to the host AND the explicit consumer stream.
// cuMemcpyHtoD from a pageable std::vector is not sufficient here: its DMA
// may still be running on the default stream when a non-blocking stream
// consumes the destination. This function never uses the default stream.
inline void upload_host_to_device_sync(const CudaContext& context,
    CUdeviceptr destination, const void* source, std::size_t bytes)
{
    if(bytes == 0) return;
    if(destination == 0 || source == nullptr)
        throw std::invalid_argument("Nonempty CUDA upload requires valid pointers.");
    detail::PinnedHostUpload staging(context, bytes);
    staging.copy_to(destination, source, bytes);
}
} // namespace rainbow
