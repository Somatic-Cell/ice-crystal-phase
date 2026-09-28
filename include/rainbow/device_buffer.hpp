#pragma once

#include <rainbow/cuda_driver.hpp>
#include <rainbow/cuda_error.hpp>

#include <cuda.h>

#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>

namespace rainbow
{


// 通常の CUDA デバイスメモリを，T の要素数とともに共有するホスト側の方

template<class T>
class DeviceBuffer final
{
    static_assert(std::is_trivially_copyable_v<T>,
    "DeviceBuffer requires a trivially-copyable element type.");
    static_assert(!std::is_const_v<T> && !std::is_volatile_v<T>,
    "DeviceBuffer owns mutable storage; T must NOT be cv-qualified.");

    // CUDA の通常の確保が保証する 256-byte alignment の範囲に限定
    static_assert(alignof(T) <= 256,
    "Over-aligned types require a different allocation strategy.");
    
public:
    // 空のバッファを作る．GPU のメモリ確保は， allocate() で明示的に行う
    explicit DeviceBuffer(const CudaContext& cuda_context) noexcept
        : cuda_context_(cuda_context)
    {
    }

    ~DeviceBuffer() noexcept
    {
        static_cast<void>(close_noexcept());
    }

    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(DeviceBuffer&) = delete;

    DeviceBuffer(DeviceBuffer&&) = delete;
    DeviceBuffer& operator=(DeviceBuffer&&) = delete;

    // 空の状態にのみ確保できる
    void allocate(const std::size_t element_count)
    {
        if(!is_empty())
        {
            throw std::logic_error("DeviceBuffer is already allocated.");
        }
        if(element_count == 0)
        {
            return;
        }
        if(element_count > (std::numeric_limits<std::size_t>::max)() / sizeof(T))
        {
            throw std::length_error("DeviceBuffer byte size would overflow.");
        }

        cuda_context_.make_current();
        CUdeviceptr new_address = 0;
        RAINBOW_CUDA_CHECK(cuMemAlloc(&new_address, element_count * sizeof(T)));

        // API が成功してから状態を更新する．以前の代入は例外を出さない
        device_address_ = new_address;
        element_count_ = element_count;
    }

    // 全領域のビットを 0 にする
    void zero_byte_async(const CUstream stream)
    {
        require_stream(stream);
        if(is_empty())
        {
            return;
        }
        cuda_context_.make_current();
        RAINBOW_CUDA_CHECK(cuMemsetD8Async(device_address_, 0, byte_size(), stream));
    }

    // ホストからバッファ全体へ転送する
    void upload_async(const std::span<const T> source, const CUstream stream)
    {
        require_element_count(source.size());
        require_stream(stream);
        if(source.empty())
        {
            return;
        }
        if(source.data() == nullptr)
        {
            throw std::invalid_argument("DeviceBuffer upload source is null.");
        }
        cuda_context_.make_current();
        RAINBOW_CUDA_CHECK(cuMemcpyHtoDAsync(
            device_address_, source.data(), byte_size(), stream));
    }

    // 全領域をホストへ同期コピーする
    void download(const std::span<T> destination) const
    {
        require_element_count(destination.size());
        if(destination.empty())
        {
            return;
        }
        if(destination.data() == nullptr)
        {
            throw std::invalid_argument("DeviceBuffer download destination is null.");
        }
        cuda_context_.make_current();
        RAINBOW_CUDA_CHECK(cuMemcpyDtoH(
            destination.data(), device_address_, byte_size()));
    }

    // 明示的に開放し，失敗したら例外で通知
    void close()
    {
        if(!close_noexcept())
        {
            throw std::runtime_error("DeviceBuffer cleanup failed. See stderr.");
        }
    }

    [[nodiscard]]
    bool close_noexcept() noexcept
    {
        if(is_empty())
        {
            return true;
        }
        if(!detail::report_cuda_cleanup_result(
               cuCtxSetCurrent(cuda_context_.handle()), "cuCtxSetCurrent(buffer cleanup)"))
        {
            return false;
        }

        const CUresult result = cuMemFree(device_address_);
        device_address_ = 0;
        element_count_ = 0;
        return detail::report_cuda_cleanup_result(result, "cuMemFree(DeviceBuffer)");
    }

    [[nodiscard]]
    bool is_empty() const noexcept { return device_address_ == 0; }

    [[nodiscard]]
    std::size_t element_count() const noexcept { return element_count_; }

    [[nodiscard]]
    std::size_t byte_size() const noexcept { return element_count_ * sizeof(T); }

    // Driver API / OptiX の host API 向けのアドレス．所有権は移さない
    [[nodiscard]]
    CUdeviceptr address() const noexcept { return device_address_; }

    // 
    [[nodiscard]]
    T* data() noexcept { return reinterpret_cast<T*>(device_address_); }

private:
    void require_element_count(const std::size_t count) const
    {
        if(count != element_count_)
        {
            throw std::invalid_argument("DeviceBuffer transfer element count mismatch.");
        }
    }

    // これは T の状態には属さない，内部専用の引数のチェック．公開の utility class は作らない
    static void require_stream(const CUstream stream)
    {
        if(stream == nullptr)
        {
            throw std::invalid_argument("DeviceBuffer requires an explicit CUDA stream.");
        }
    }

    const CudaContext& cuda_context_;
    CUdeviceptr device_address_ = 0;
    std::size_t element_count_ = 0;

};

} // namespace rainbow