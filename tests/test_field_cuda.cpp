#include "field_test_cases.hpp"

#include <rainbow/cuda_driver.hpp>
#include <rainbow/cuda_error.hpp>
#include <rainbow/cuda_module.hpp>
#include <rainbow/device_buffer.hpp>

#include <cuda.h>

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace
{

// 所有する状態と，非同期処理の完了待ちに責任を持つテスト実行オブジェクト．
// メモリ確保・module load・ファイル読込み・エラー処理の本体は既存クラスを使う．
class FieldCudaTest final
{
public:
    explicit FieldCudaTest(const rainbow::CudaContext& cuda_context)
        : cuda_context_(cuda_context),
          inputs_(rainbow::tests::make_field_test_inputs()),
          outputs_(inputs_.size()),
          module_(cuda_context),
          device_inputs_(cuda_context),
          device_outputs_(cuda_context)
    {
    }

    ~FieldCudaTest() noexcept
    {
        // destructor 本体はメンバの破棄より前．転送元 inputs_ が生存中に同期する．
        static_cast<void>(release_resources());
    }

    FieldCudaTest(const FieldCudaTest&) = delete;
    FieldCudaTest& operator=(const FieldCudaTest&) = delete;
    FieldCudaTest(FieldCudaTest&&) = delete;
    FieldCudaTest& operator=(FieldCudaTest&&) = delete;

    void run(const std::filesystem::path& fatbin_path)
    {
        if(inputs_.size() > (std::numeric_limits<std::uint32_t>::max)())
        {
            throw std::length_error("Too many field test cases.");
        }
        cuda_context_.make_current();
        module_.load_fatbin(fatbin_path);
        const CUfunction function = module_.find_function("evaluate_field_arithmetic");
        device_inputs_.allocate(inputs_.size());
        device_outputs_.allocate(outputs_.size());

        const CUstream stream = cuda_context_.stream();
        // API が途中で失敗しても destructor で完了待ちを試みるため，先に立てる．
        has_pending_work_ = true;
        device_inputs_.upload_async(
            std::span<const rainbow::tests::FieldTestInput>{inputs_}, stream);
        // 最新リポジトリのメンバ名 (zero_byte_async) に合わせる．
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

        device_outputs_.download(std::span<rainbow::tests::FieldTestResult>{outputs_});
        rainbow::tests::FieldTestVerifier verifier;
        verifier.verify(
            std::span<const rainbow::tests::FieldTestInput>{inputs_},
            std::span<const rainbow::tests::FieldTestResult>{outputs_});
        if(!release_resources())
        {
            throw std::runtime_error("Field CUDA test cleanup failed. See stderr.");
        }
        verifier.print_summary("CUDA", inputs_.size());
    }

private:
    [[nodiscard]]
    bool release_resources() noexcept
    {
        bool succeeded = true;
        if(has_pending_work_)
        {
            if(!rainbow::detail::report_cuda_cleanup_result(
                   cuCtxSetCurrent(cuda_context_.handle()), "cuCtxSetCurrent(field test)"))
            {
                return false;
            }
            succeeded = rainbow::detail::report_cuda_cleanup_result(
                cuStreamSynchronize(cuda_context_.stream()), "cuStreamSynchronize(field test)");
            has_pending_work_ = false;
        }
        // API 側を && の左に置き，一つが失敗しても後続の解放を省略しない．
        succeeded = device_outputs_.close_noexcept() && succeeded;
        succeeded = device_inputs_.close_noexcept() && succeeded;
        succeeded = module_.close_noexcept() && succeeded;
        return succeeded;
    }

    // context を借用する．main() の owner がこのオブジェクトより長く生存する．
    const rainbow::CudaContext& cuda_context_;
    std::vector<rainbow::tests::FieldTestInput> inputs_;
    std::vector<rainbow::tests::FieldTestResult> outputs_;
    rainbow::CudaModule module_;
    rainbow::DeviceBuffer<rainbow::tests::FieldTestInput> device_inputs_;
    rainbow::DeviceBuffer<rainbow::tests::FieldTestResult> device_outputs_;
    bool has_pending_work_ = false;
};

// コマンドライン文字型の違いを入口だけに閉じ込める．GPU ロジックは共通．
int run_test(const std::filesystem::path& fatbin_path)
{
    try
    {
        rainbow::CudaContext cuda_context{0};
        FieldCudaTest test(cuda_context);
        test.run(fatbin_path);
        return EXIT_SUCCESS;
    }
    catch(const std::exception& exception)
    {
        std::cerr << "Error: " << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}

} // namespace

// Windows では CTest が渡す日本語等のパスもワイド文字として受け取る．
// 既存 apps/main.cpp は変更しない．このファイルは別のテスト executable の入口．
#if defined(_WIN32)
int wmain(const int argc, wchar_t* argv[])
#else
int main(const int argc, char* argv[])
#endif
{
    if(argc != 2)
    {
        std::cerr << "Usage: rainbow_field_cuda_tests <field_test.fatbin>\n";
        return EXIT_FAILURE;
    }
    try
    {
        return run_test(std::filesystem::path{argv[1]});
    }
    catch(const std::exception& exception)
    {
        std::cerr << "Error constructing module path: " << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
