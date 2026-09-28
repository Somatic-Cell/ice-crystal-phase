#pragma once

#include <rainbow/cuda_driver.hpp>
#include <rainbow/cuda_error.hpp>
#include <rainbow/cuda_module.hpp>
#include <rainbow/device_buffer.hpp>

#include <cuda.h>

#include <cstdint>
#include <filesystem>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace rainbow::tests
{

// field と phase の両テストに実在する共通処理だけを集約する．本番 kernel の
// 汎用フレームワークではなく，(const Input*, Output*, uint32_t) という検査用
// kernel ABI に限定する．期待値生成・比較はこのクラスの責務ではない．
//
// 入力を所有する理由: 非同期転送後に例外が起きても，destructor の完了待ちまで
// ホスト転送元を生存させるため．context は借用し，ここでは破棄しない．
template<class Input, class Output>
class CudaArrayTest final
{
public:
    explicit CudaArrayTest(const CudaContext& cuda_context, std::vector<Input> inputs)
        : cuda_context_(cuda_context), inputs_(std::move(inputs)), outputs_(inputs_.size()),
          module_(cuda_context), device_inputs_(cuda_context), device_outputs_(cuda_context)
    {
    }

    ~CudaArrayTest() noexcept
    {
        static_cast<void>(release_resources());
    }

    CudaArrayTest(const CudaArrayTest&) = delete;
    CudaArrayTest& operator=(const CudaArrayTest&) = delete;
    CudaArrayTest(CudaArrayTest&&) = delete;
    CudaArrayTest& operator=(CudaArrayTest&&) = delete;

    void run(const std::filesystem::path& fatbin_path, const char* kernel_name)
    {
        if(inputs_.empty() || inputs_.size() > (std::numeric_limits<std::uint32_t>::max)())
        {
            throw std::length_error("CUDA array test requires 1..UINT32_MAX cases.");
        }
        cuda_context_.make_current();
        module_.load_fatbin(fatbin_path);
        const CUfunction function = module_.find_function(kernel_name);
        device_inputs_.allocate(inputs_.size());
        device_outputs_.allocate(outputs_.size());

        const CUstream stream = cuda_context_.stream();
        has_pending_work_ = true;
        device_inputs_.upload_async(std::span<const Input>{inputs_}, stream);
        device_outputs_.zero_byte_async(stream);

        CUdeviceptr input_address = device_inputs_.address();
        CUdeviceptr output_address = device_outputs_.address();
        std::uint32_t count = static_cast<std::uint32_t>(inputs_.size());
        void* arguments[] = {&input_address, &output_address, &count};
        constexpr unsigned int threads_per_block = 256;
        const unsigned int block_count = count / threads_per_block
            + (count % threads_per_block != 0 ? 1u : 0u);
        RAINBOW_CUDA_CHECK(cuLaunchKernel(
            function, block_count, 1, 1, threads_per_block, 1, 1,
            0, stream, arguments, nullptr));
        RAINBOW_CUDA_CHECK(cuStreamSynchronize(stream));
        has_pending_work_ = false;
        device_outputs_.download(std::span<Output>{outputs_});

        // 結果の検証より先に GPU 資源の解放を完了する．ホスト出力は保持される．
        if(!release_resources())
        {
            throw std::runtime_error("CUDA array test cleanup failed. See stderr.");
        }
    }

    [[nodiscard]]
    std::span<const Input> inputs() const noexcept { return inputs_; }

    // 正常に run() が戻った後で読む．返す span の所有者はこのオブジェクト．
    [[nodiscard]]
    std::span<const Output> outputs() const noexcept { return outputs_; }

private:
    [[nodiscard]]
    bool release_resources() noexcept
    {
        bool succeeded = true;
        if(has_pending_work_)
        {
            if(!detail::report_cuda_cleanup_result(
                   cuCtxSetCurrent(cuda_context_.handle()), "cuCtxSetCurrent(array test)"))
            {
                return false;
            }
            succeeded = detail::report_cuda_cleanup_result(
                cuStreamSynchronize(cuda_context_.stream()), "cuStreamSynchronize(array test)");
            has_pending_work_ = false;
        }
        // 同期後にバッファ・module を破棄する．context 喪失からの復旧は保証しない．
        succeeded = device_outputs_.close_noexcept() && succeeded;
        succeeded = device_inputs_.close_noexcept() && succeeded;
        succeeded = module_.close_noexcept() && succeeded;
        return succeeded;
    }

    const CudaContext& cuda_context_;
    std::vector<Input> inputs_;
    std::vector<Output> outputs_;
    CudaModule module_;
    DeviceBuffer<Input> device_inputs_;
    DeviceBuffer<Output> device_outputs_;
    bool has_pending_work_ = false;
};

} // namespace rainbow::tests
