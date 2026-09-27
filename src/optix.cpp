#include <rainbow/optix.hpp>

#include <rainbow/cuda_driver.hpp>
#include <rainbow/trace_launch_params.hpp>

#include <cuda.h>
#include <optix.h>

#include <optix_function_table_definition.h>
#include <optix_stack_size.h>
#include <optix_stubs.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace
{

// エラー処理
// 通常処理: 戻り値を検査し， API名，ファイル，行番号を付けた例外を出す
// 廃棄処理: 元の例外を上書きしないため， stderr への記録にとどめる
void check_cuda(
    const CUresult result,
    const char* expression,
    const char* file,
    const int line)
{
    if(result == CUDA_SUCCESS)
    {
        return;
    }

    const char* error_name = nullptr;
    const char* error_description = nullptr;
    static_cast<void>(cuGetErrorName(result, &error_name));
    static_cast<void>(cuGetErrorString(result, &error_description));

    std::ostringstream message;

    message 
        << expression
        << " failed at "
        << file 
        << " : "
        << line
        << "\n"
        << "CUDA Error: "
        <<  (error_name != nullptr ? error_name : "unknown")
        << " (" << static_cast<int>(result) << ")"
        << "\n"
        << "Description: "
        << (error_description != nullptr ? error_description : "unavailable");
    throw std::runtime_error(message.str());
}

void check_optix(
    const OptixResult result,
    const char* expression,
    const char* file,
    const int line)
{
    if(result == OPTIX_SUCCESS)
    {
        return;
    }

    // optix_stubs.h は，optixInit() の失敗時にもエラー名を取得できる実装になっている
    const char* error_name = optixGetErrorName(result);
    const char* error_description = optixGetErrorString(result);
    
    std::ostringstream message;
    message
        << expression
        << " failed at "
        << file
        << ':'
        << line
        << '\n'
        << "OptiX error: "
        << (error_name != nullptr ? error_name : "unknown")
        << " (" << static_cast<int>(result) << ")"
        << '\n'
        << "Description: "
        << (error_description != nullptr ? error_description : "unavailable");
    throw std::runtime_error(message.str());
}

#define RAINBOW_M2_CUDA_CHECK(expression) \
    check_cuda((expression), #expression, __FILE__, __LINE__)

#define RAINBOW_M2_OPTIX_CHECK(expression)                         \
    check_optix((expression), #expression, __FILE__, __LINE__)


[[nodiscard]]
bool report_cleanup_result(
    const CUresult result,
    const char* expression) noexcept
{
    if(result == CUDA_SUCCESS)
    {
        return true;
    }

    const char* error_name = nullptr;
    static_cast<void>(cuGetErrorName(result, &error_name));
    std::fprintf(
        stderr, "[CUDA cleanup] %s: %s: (%d)\n",
        expression,
        error_name!= nullptr ? error_name : "unknown",
        static_cast<int>(result)
    );
    return false;
}

[[nodiscard]]
bool report_cleanup_result(
    const OptixResult result,
    const char* expression) noexcept
{
    if(result == OPTIX_SUCCESS)
    {
        return true;
    }

    const char* error_name = optixGetErrorName(result);
    std::fprintf(
        stderr, "[OptiX cleanup] %s: %s: (%d)\n",
        expression,
        error_name != nullptr ? error_name : "unknown",
        static_cast<int>(result)
    );
    return false;
}

// -------------------------
// バイナリファイルを std::vector<char> に読み込んで返す
// -------------------------
[[nodiscard]]
std::vector<char> read_binary_file(const std::filesystem::path& path)
{
    // ファイルの末尾を調べたいため
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if(!input)
    {
        throw std::runtime_error("Failed to open Optix IR: " + path.string());
    }

    const std::streampos end_position = input.tellg();
    if(end_position == std::streampos(-1))
    {
        throw std::runtime_error("Failed to read OptiX IR size: " + path.string());
    }
    
    // 空のファイルは拒否
    const std::streamoff file_byte_size =
        static_cast<std::streamoff>(end_position);
    if(file_byte_size <= 0)
    {
        throw std::runtime_error("OptiX IR is empty: " + path.string());
    }

    // std::uintmax_t で表現しきれないサイズも拒否
    const std::uintmax_t byte_size = static_cast<std::uintmax_t>(file_byte_size);
    if(
        byte_size > static_cast<std::uintmax_t>(std::numeric_limits<std::size_t>::max())
        || byte_size > static_cast<std::uintmax_t>(std::numeric_limits<std::streamsize>::max())
    )
    {
        throw std::runtime_error("Optix IR is too large: " + path.string());
    }

    // 先頭から byte_size 分のデータを読む
    std::vector<char> bytes(static_cast<std::size_t>(byte_size));
    input.seekg(0, std::ios::beg);
    if(!input.read(bytes.data(), static_cast<std::streamsize>(bytes.size())))
    {
        throw std::runtime_error("Failed to read OptiX IR: " + path.string());
    }

    return bytes;
}


void optix_log_callback(
    const unsigned int level,
    const char* tag,
    const char* message,
    void*) noexcept
{
    // SDK の callback から C++ へ例外を出さない
    std::fprintf(
        stderr, 
        "[OptiX][%u][%s] %s\n",
        level,
        tag != nullptr ? tag : "",
        message != nullptr ? message : ""
    );
}

constexpr std::size_t log_buffer_byte_size = 16384;
using OptixLogBuffer = std::array<char, log_buffer_byte_size>;

void print_optix_log(
    const char* stage,
    const OptixLogBuffer& log,
    const std::size_t reported_byte_size) noexcept
{
    const auto end = std::find(log.begin(), log.end(), '\0');
    const auto stored_byte_size = static_cast<std::size_t>(end - log.begin());
    if(stored_byte_size != 0)
    {
        std::fprintf(stderr, "[OptiX %s log]\n", stage);
        static_cast<void>(std::fwrite(log.data(), 1, stored_byte_size, stderr));
        std::fputc('\n', stderr);
    }
    if(reported_byte_size > log.size())
    {
        std::fprintf(
            stderr,
            "[OptiX %s log truncated: requrired %zu bytes, capacity %zu bytes]\n",
            stage, reported_byte_size, log.size()
        );
    }
}

void initialize_optix()
{
    static std::once_flag initialization_flag;
    std::call_once(initialization_flag, []
    {
        RAINBOW_M2_OPTIX_CHECK(optixInit());
    });
}

// -----------------------
// SBT / launch parameters のレイアウト
// -----------------------

struct alignas(OPTIX_SBT_RECORD_ALIGNMENT) RaygenRecord
{
    char header[OPTIX_SBT_RECORD_HEADER_SIZE];
};

struct alignas(OPTIX_SBT_RECORD_ALIGNMENT) MissRecord
{
    char header[OPTIX_SBT_RECORD_HEADER_SIZE];
};

static_assert(sizeof(void*) == 8, "This M2 implementation requires a 64-bit build.");
static_assert(sizeof(CUdeviceptr) == sizeof(void*));
static_assert(sizeof(unsigned int) == sizeof(std::uint32_t));
static_assert(alignof(RaygenRecord) == OPTIX_SBT_RECORD_ALIGNMENT);
static_assert(sizeof(RaygenRecord) % OPTIX_SBT_RECORD_ALIGNMENT == 0);
static_assert(std::is_trivially_copyable_v<RaygenRecord>);
static_assert(std::is_standard_layout_v<rainbow::TraceLaunchParams>);
static_assert(std::is_trivially_copyable_v<rainbow::TraceLaunchParams>);
static_assert(sizeof(rainbow::TraceLaunchParams) == 16);
static_assert(alignof(rainbow::TraceLaunchParams) == 8);
static_assert(offsetof(rainbow::TraceLaunchParams, output) == 0);
static_assert(offsetof(rainbow::TraceLaunchParams, count) == 8);
static_assert(offsetof(rainbow::TraceLaunchParams, reserved) == 12);

class OptixSmokeTest final
{
public:
    explicit OptixSmokeTest(const rainbow::CudaContext& cuda_context) noexcept
        : cuda_context_(cuda_context)
    {
    }

    ~OptixSmokeTest() noexcept
    {
        // 例外経路の安全網．失敗はログへ記録し，元の例外を保持する．
        static_cast<void>(release_resources());
    }

    OptixSmokeTest(const OptixSmokeTest&) = delete;
    OptixSmokeTest& operator=(const OptixSmokeTest&) = delete;
    OptixSmokeTest(OptixSmokeTest&&) = delete;
    OptixSmokeTest& operator=(OptixSmokeTest&&) = delete;

    // この一時オブジェクトに対して一回だけ呼ぶ．
    void run(const std::filesystem::path& optixir_path);

private:
    [[nodiscard]]
    bool release_resources() noexcept;

    const rainbow::CudaContext& cuda_context_;

    OptixDeviceContext optix_context_ = nullptr;
    OptixModule module_ = nullptr;
    OptixProgramGroup raygen_program_group_ = nullptr;
    OptixProgramGroup miss_program_group_ = nullptr;
    OptixPipeline pipeline_ = nullptr;

    CUdeviceptr device_output_ = 0;
    CUdeviceptr device_params_ = 0;
    CUdeviceptr device_raygen_record_ = 0;
    CUdeviceptr device_miss_record_ = 0;

    // Async 転送元をローカルな一時変数にせず，所有者のメンバとして保持する．
    // 例外が起きても destructor 本体の同期が終わるまで，転送元は生存している．
    rainbow::TraceLaunchParams host_params_ = {};
    RaygenRecord host_raygen_record_ = {};
    MissRecord host_miss_record_ = {};
    bool has_pending_work_ = false;
};

void OptixSmokeTest::run(const std::filesystem::path& optixir_path)
{
    cuda_context_.make_current();
    if(cuda_context_.handle() == nullptr || cuda_context_.stream() == nullptr)
    {
        throw std::invalid_argument(
            "M2 requires a valid CudaContext and an explicitly created CUDA stream.");
    }

    const std::vector<char> optixir = read_binary_file(optixir_path);
    initialize_optix();

    // ----- OptiX device context: CUDA context は新規作成しない． -----
    OptixDeviceContextOptions context_options = {};
    context_options.logCallbackFunction = &optix_log_callback;
    context_options.logCallbackLevel = 4;
#ifndef NDEBUG
    context_options.validationMode = OPTIX_DEVICE_CONTEXT_VALIDATION_MODE_ALL;
#endif
    RAINBOW_M2_OPTIX_CHECK(optixDeviceContextCreate(
        cuda_context_.handle(), &context_options, &optix_context_));

    // ----- IR -> module -----
    OptixModuleCompileOptions module_options = {};
    module_options.maxRegisterCount = OPTIX_COMPILE_DEFAULT_MAX_REGISTER_COUNT;
    module_options.optLevel = OPTIX_COMPILE_OPTIMIZATION_DEFAULT;
    module_options.debugLevel = OPTIX_COMPILE_DEBUG_LEVEL_DEFAULT;

    OptixPipelineCompileOptions pipeline_compile_options = {};
    pipeline_compile_options.usesMotionBlur = 0;
    pipeline_compile_options.traversableGraphFlags =
        OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_GAS;
    pipeline_compile_options.numPayloadValues = 0;
    pipeline_compile_options.numAttributeValues = 0;
    pipeline_compile_options.exceptionFlags = OPTIX_EXCEPTION_FLAG_NONE;
    pipeline_compile_options.pipelineLaunchParamsVariableName = "params";
    pipeline_compile_options.pipelineLaunchParamsSizeInBytes =
    sizeof(rainbow::TraceLaunchParams);
    pipeline_compile_options.usesPrimitiveTypeFlags = OPTIX_PRIMITIVE_TYPE_FLAGS_CUSTOM;
    // 上の graph / primitive flags はコンパイル上の設定で，GAS を作る処理ではない．

    OptixLogBuffer log = {};
    std::size_t log_byte_size = log.size();
    OptixResult result = optixModuleCreate(
        optix_context_, &module_options, &pipeline_compile_options,
        optixir.data(), optixir.size(), log.data(), &log_byte_size, &module_);
    print_optix_log("module", log, log_byte_size);
    check_optix(result, "optixModuleCreate", __FILE__, __LINE__);

    // ----- module 内の raygen を program group として指定する． -----
    OptixProgramGroupDesc raygen_description = {};
    raygen_description.kind = OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
    raygen_description.raygen.module = module_;
    raygen_description.raygen.entryFunctionName = "__raygen__raindrop_trace";

    OptixProgramGroupDesc miss_description = {};
    miss_description.kind = OPTIX_PROGRAM_GROUP_KIND_MISS;

    // M2 では実行する miss shader を指定しない．
    miss_description.miss.module = nullptr;
    miss_description.miss.entryFunctionName = nullptr;

    OptixProgramGroupOptions program_group_options = {};
    log.fill('\0');
    log_byte_size = log.size();
    result = optixProgramGroupCreate(
        optix_context_, &raygen_description, 1, &program_group_options,
        log.data(), &log_byte_size, &raygen_program_group_);
    print_optix_log("program group", log, log_byte_size);
    check_optix(result, "optixProgramGroupCreate", __FILE__, __LINE__);

    // ----- raygen-only pipeline -----
    OptixPipelineLinkOptions pipeline_link_options = {};
    // M2 は optixTrace() を一度も呼ばないので，許容する trace 深度は 0 にする．
    // 雨粒追跡を追加する時点で，その呼び出し構造に合わせて変更する．
    pipeline_link_options.maxTraceDepth = 0;

    log.fill('\0');
    log_byte_size = log.size();
    result = optixPipelineCreate(
        optix_context_, &pipeline_compile_options, &pipeline_link_options,
        &raygen_program_group_, 1, log.data(), &log_byte_size, &pipeline_);
    print_optix_log("pipeline", log, log_byte_size);
    check_optix(result, "optixPipelineCreate", __FILE__, __LINE__);

    // ----- 必要 stack size を program group の情報から求める． -----
    OptixStackSizes stack_sizes = {};
    RAINBOW_M2_OPTIX_CHECK(optixUtilAccumulateStackSizes(
        raygen_program_group_, &stack_sizes, pipeline_));

    unsigned int direct_callable_from_traversal = 0;
    unsigned int direct_callable_from_state = 0;
    unsigned int continuation_stack_byte_size = 0;
    RAINBOW_M2_OPTIX_CHECK(optixUtilComputeStackSizes(
        &stack_sizes, pipeline_link_options.maxTraceDepth, 0, 0,
        &direct_callable_from_traversal, &direct_callable_from_state,
        &continuation_stack_byte_size));

    // 最後の 1 は traversable graph depth の上限で，trace 再帰深度とは別物．
    // M2 では traversable を使わないが，設定上は single GAS の上限に合わせる．
    RAINBOW_M2_OPTIX_CHECK(optixPipelineSetStackSize(
        pipeline_, direct_callable_from_traversal, direct_callable_from_state,
        continuation_stack_byte_size, 1));

    // ----- GPU buffer と raygen SBT record -----
    constexpr std::uint32_t element_count = 1024;
    constexpr std::uint32_t expected_mask = 0xa5a5a5a5u;
    const std::size_t output_byte_size =
        static_cast<std::size_t>(element_count) * sizeof(std::uint32_t);

    RAINBOW_M2_CUDA_CHECK(cuMemAlloc(&device_output_, output_byte_size));
    RAINBOW_M2_CUDA_CHECK(cuMemAlloc(&device_params_, sizeof(host_params_)));
    RAINBOW_M2_CUDA_CHECK(cuMemAlloc(
        &device_raygen_record_, sizeof(host_raygen_record_)));

    host_params_.output = reinterpret_cast<std::uint32_t*>(device_output_);
    host_params_.count = element_count;
    host_params_.reserved = 0;
    RAINBOW_M2_OPTIX_CHECK(optixSbtRecordPackHeader(
        raygen_program_group_, &host_raygen_record_));

    // ----- 同じ stream に初期化・転送・launch を順に投入する． -----
    // CUDA 13.x では pageable host memory を使う Async 転送も扱えるが，
    // 内部 staging によりホストをブロックする場合がある．M2 は転送速度を測らない．
    // 必要になれば pinned staging buffer へ変更する．現在も転送元の寿命は保持する．
    const CUstream stream = cuda_context_.stream();
    has_pending_work_ = true;
    RAINBOW_M2_CUDA_CHECK(cuMemsetD32Async(
        device_output_, 0, element_count, stream));
    RAINBOW_M2_CUDA_CHECK(cuMemcpyHtoDAsync(
        device_params_, &host_params_, sizeof(host_params_), stream));
    RAINBOW_M2_CUDA_CHECK(cuMemcpyHtoDAsync(
        device_raygen_record_, &host_raygen_record_, sizeof(host_raygen_record_), stream));

    OptixShaderBindingTable sbt = {};
    sbt.raygenRecord = device_raygen_record_;
    // optixTrace() がないため，miss / hitgroup / callable record は不要．
    RAINBOW_M2_OPTIX_CHECK(optixLaunch(
        pipeline_, stream, device_params_, sizeof(host_params_), &sbt,
        element_count, 1, 1));

    // launch 自体の成功と，GPU 上での実行完了は別なので，ここでも戻り値を検査する．
    RAINBOW_M2_CUDA_CHECK(cuStreamSynchronize(stream));
    has_pending_work_ = false;

    // GPU 書き込み完了後に，同期コピーで通常のホストメモリへ読み戻す．
    std::vector<std::uint32_t> host_output(element_count);
    RAINBOW_M2_CUDA_CHECK(cuMemcpyDtoH(
        host_output.data(), device_output_, output_byte_size));

    for(std::uint32_t element_index = 0; element_index < element_count; ++element_index)
    {
        const std::uint32_t expected = element_index ^ expected_mask;
        if(host_output[element_index] != expected)
        {
            std::ostringstream message;
            message << "OptiX output mismatch at index " << element_index
                    << ": expected " << expected
                    << ", got " << host_output[element_index]
                    << "\nCheck the shader ABI, raygen entry, and the staged .optixir file.";
            throw std::runtime_error(message.str());
        }
    }

    // 正常終了時の cleanup 失敗はログだけにせず，呼び出し元へ通知する．
    // 例外中の destructor では投げない，という経路との区別を明確にする．
    if(!release_resources())
    {
        throw std::runtime_error(
            "OptiX output passed verification, but resource cleanup failed. See stderr.");
    }
    std::cout << "OptiX raygen smoke test: passed\n";
}

bool OptixSmokeTest::release_resources() noexcept
{
    if(optix_context_ == nullptr && module_ == nullptr
       && raygen_program_group_ == nullptr && pipeline_ == nullptr
       && device_output_ == 0 && device_params_ == 0 && device_raygen_record_ == 0)
    {
        return true;
    }

    // 借りた context を current にするが，所有者の stream / retain 参照は解放しない．
    if(!report_cleanup_result(
           cuCtxSetCurrent(cuda_context_.handle()), "cuCtxSetCurrent"))
    {
        // 異なる context のままハンドルを解放してはいけない．
        // driver/context 自体が失われた場合，この destructor だけでは復旧できない．
        return false;
    }

    bool succeeded = true;
    if(has_pending_work_)
    {
        succeeded = report_cleanup_result(
            cuStreamSynchronize(cuda_context_.stream()), "cuStreamSynchronize") && succeeded;
        has_pending_work_ = false;
    }

    // 短絡評価で後続 cleanup を飛ばさないよう，API 呼び出しを && の左側に置く．
    if(device_raygen_record_ != 0)
    {
        succeeded = report_cleanup_result(
            cuMemFree(device_raygen_record_), "cuMemFree(raygen record)") && succeeded;
        device_raygen_record_ = 0;
    }
    if(device_params_ != 0)
    {
        succeeded = report_cleanup_result(
            cuMemFree(device_params_), "cuMemFree(launch params)") && succeeded;
        device_params_ = 0;
    }
    if(device_output_ != 0)
    {
        succeeded = report_cleanup_result(
            cuMemFree(device_output_), "cuMemFree(output)") && succeeded;
        device_output_ = 0;
    }
    if(pipeline_ != nullptr)
    {
        succeeded = report_cleanup_result(
            optixPipelineDestroy(pipeline_), "optixPipelineDestroy") && succeeded;
        pipeline_ = nullptr;
    }
    if(raygen_program_group_ != nullptr)
    {
        succeeded = report_cleanup_result(
            optixProgramGroupDestroy(raygen_program_group_), "optixProgramGroupDestroy")
            && succeeded;
        raygen_program_group_ = nullptr;
    }
    if(module_ != nullptr)
    {
        succeeded = report_cleanup_result(
            optixModuleDestroy(module_), "optixModuleDestroy") && succeeded;
        module_ = nullptr;
    }
    if(optix_context_ != nullptr)
    {
        succeeded = report_cleanup_result(
            optixDeviceContextDestroy(optix_context_), "optixDeviceContextDestroy") && succeeded;
        optix_context_ = nullptr;
    }
    return succeeded;
}

} // namespace

namespace rainbow
{

void run_optix_smoke_test(
    const CudaContext& cuda_context,
    const std::filesystem::path& optixir_path)
{
    // cuda_context は呼び出し元が所有し，この一時オブジェクトより長く生存する．
    // 外部 CUDA context の切替えを復元する汎用 API ではなく，M2 の単一スレッド用．
    OptixSmokeTest smoke_test(cuda_context);
    smoke_test.run(optixir_path);
}

} // namespace rainbow