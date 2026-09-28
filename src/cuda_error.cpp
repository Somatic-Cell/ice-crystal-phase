#include <rainbow/cuda_error.hpp>

#include <cstdio>
#include <sstream>
#include <stdexcept>

namespace rainbow::detail
{
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

    // エラーの取得に失敗しても，元のエラーコードによる診断は残す
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

bool report_cuda_cleanup_result(
    const CUresult result,
    const char* expression) noexcept
{
    if(result == CUDA_SUCCESS)
    {
        return true;
    }

    // デストラクタから使う想定のため，例外を投げる関数はよばない
    // ここで扱うのは解法処理の結果で，リソースの解法自体は行わない
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


} // namespace rainbow::detail