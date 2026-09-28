#include <rainbow/cuda_driver.hpp>
#include <rainbow/cuda_error.hpp>
#include <rainbow/cuda_module.hpp>
#include <rainbow/device_buffer.hpp>

#include <cuda.h>

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <span>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace 
{


class CudaDriverSmokeTest final
{
public:
    explicit CudaDriverSmokeTest(const rainbow::CudaContext& cuda_context) noexcept
        : cuda_context_(cuda_context),
        module_(cuda_context),
        device_output_(cuda_context)
    {
    }

    ~CudaDriverSmokeTest() noexcept
    {
        // 通常経路で解放済みなら何もしないが，例外経路では cleanup の安全網になる
        static_cast<void>(release_resources());
    }

    CudaDriverSmokeTest(const CudaDriverSmokeTest&) = delete;
    CudaDriverSmokeTest& operator=(const CudaDriverSmokeTest&) = delete;
    CudaDriverSmokeTest(CudaDriverSmokeTest&&) = delete;
    CudaDriverSmokeTest& operator=(CudaDriverSmokeTest&&) = delete;

    void run(const std::filesystem::path& fatbin_path);

private:
    [[nodiscard]]
    bool release_resources() noexcept;

    const rainbow::CudaContext& cuda_context_;
    rainbow::CudaModule module_;
    rainbow::DeviceBuffer<std::uint32_t> device_output_;
    bool has_pending_work_ = false;
};

void CudaDriverSmokeTest::run(const std::filesystem::path& fatbin_path)
{
    cuda_context_.make_current();
    module_.load_fatbin(fatbin_path);
    const CUfunction function = module_.find_function("write_test_pattern");

    std::uint32_t element_count = 1024;
    constexpr std::uint32_t expected_mask = 0x5a5a5a5au;
    constexpr unsigned int threads_per_block = 256;
    const unsigned int block_count =
        element_count / threads_per_block
        + (element_count % threads_per_block != 0 ? 1u : 0u);

    static_assert(sizeof(void*) == 8, "This project requires a 64-bit build.");
    static_assert(sizeof(CUdeviceptr) == sizeof(void*));
    static_assert(sizeof(unsigned int) == sizeof(std::uint32_t));

    device_output_.allocate(element_count);
    const CUstream stream = cuda_context_.stream();
    has_pending_work_ = true;
    device_output_.zero_byte_async(stream);

    CUdeviceptr output_address = device_output_.address();
    void* kernel_parameters[] = {&output_address, &element_count};
    RAINBOW_CUDA_CHECK(cuLaunchKernel(
        function,
        block_count, 1, 1,
        threads_per_block, 1, 1,
        0, stream, kernel_parameters, nullptr));
    
    RAINBOW_CUDA_CHECK(cuStreamSynchronize(stream));
    has_pending_work_ = false;

    std::vector<std::uint32_t> host_output(element_count);
    device_output_.download(std::span<std::uint32_t>{host_output});
    for(std::uint32_t element_index = 0; element_index < element_count; ++element_index)
    {
        const std::uint32_t expected = element_index ^ expected_mask;
        if(host_output[element_index] != expected)
        {
            std::ostringstream message;
            message << "CUDA output mismatch at index " << element_index
                    << ": expected " << expected
                    << ", got " << host_output[element_index]
                    << "\nCheck the kernel body and the staged .fatbin file.";
            throw std::runtime_error(message.str());
        }
    }

    char device_name[256] = {};
    RAINBOW_CUDA_CHECK(cuDeviceGetName(
        device_name, static_cast<int>(sizeof(device_name)), cuda_context_.device()));

    if(!release_resources())
    {
        throw std::runtime_error(
            "CUDA output passed verification, but resource cleanup failed. See stderr.");
    }
    std::cout << "CUDA device: " << device_name << '\n'
              << "CUDA Driver API smoke test: passed\n";

}

bool CudaDriverSmokeTest::release_resources() noexcept
{
    if(!module_.is_loaded() && device_output_.is_empty())
    {
        return true;
    }

    if(!rainbow::detail::report_cuda_cleanup_result(
           cuCtxSetCurrent(cuda_context_.handle()), "cuCtxSetCurrent(smoke cleanup)"))
    {
        // 間違った context のまま解放しない．driver/context 喪失時の復旧は対象外．
        return false;
    }

    bool succeeded = true;
    if(has_pending_work_)
    {
        succeeded = rainbow::detail::report_cuda_cleanup_result(
            cuStreamSynchronize(cuda_context_.stream()), "cuStreamSynchronize(smoke cleanup)")
            && succeeded;
        has_pending_work_ = false;
    }

    succeeded = device_output_.close_noexcept() && succeeded;
    succeeded = module_.close_noexcept() && succeeded;
    return succeeded;
}
}


namespace rainbow
{
CudaContext::CudaContext(
    const int device_ordinal
)
{
    // 初期化
    RAINBOW_CUDA_CHECK(
        cuInit(0)
    );

    // 入力値が不正でないかどうか確認
    int device_count = 0;

    RAINBOW_CUDA_CHECK(
        cuDeviceGetCount(&device_count)
    );

    if(device_ordinal < 0 || device_ordinal >= device_count)
    {
        throw std::runtime_error(
            "Invarid CUDA device ordinal."
        );
    }

    RAINBOW_CUDA_CHECK(
        cuDeviceGet(
            &device_,
            device_ordinal
        )
    );

    RAINBOW_CUDA_CHECK(
        cuDevicePrimaryCtxRetain(
            &context_,
            device_
        )
    );

    primary_context_retained_ = true;

    // リソースを順番に確保する．失敗したら開放し，各メンバ変数を再初期化する．
    try{
        RAINBOW_CUDA_CHECK(cuCtxSetCurrent(context_));
        RAINBOW_CUDA_CHECK(cuStreamCreate(&stream_, CU_STREAM_NON_BLOCKING));
    }
    catch(...)
    {
        if(context_ != nullptr)
        {
            static_cast<void>(
                cuCtxSetCurrent(nullptr)
            );
        }
        if(primary_context_retained_)
        {
            static_cast<void>(
                cuDevicePrimaryCtxRelease(device_)
            );
        }

        context_ = nullptr;
        primary_context_retained_ = false;

        throw;
    }
}

CudaContext::~CudaContext() noexcept
{
    if(context_ != nullptr)
    {
        static_cast<void>(cuCtxSetCurrent(context_)); // [[nodiscard]] の警告を出させないで設定
    }

    if(stream_ != nullptr)
    {
        static_cast<void>(cuStreamSynchronize(stream_));
        static_cast<void>(cuStreamDestroy(stream_));
    }

    if(context_ != nullptr)
    {
        static_cast<void>(cuCtxSetCurrent(nullptr));
    }

    if(primary_context_retained_)
    {
        static_cast<void>(cuDevicePrimaryCtxRelease(device_));
    }
}

void CudaContext::make_current() const
{
    RAINBOW_CUDA_CHECK(
        cuCtxSetCurrent(context_)
    );
}


void run_cuda_driver_smoke_test(
    const CudaContext& cuda_context,
    const std::filesystem::path& fatbin_path)
{
    CudaDriverSmokeTest smoke_test(cuda_context);
    smoke_test.run(fatbin_path);
}

} // namespace rainbow

#undef RAINBOW_CUDA_CHECK