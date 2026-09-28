#pragma once 

#include <cuda.h>

namespace rainbow::detail
{

// エラー処理
// 通常処理: 戻り値を検査し， API名，ファイル，行番号を付けた例外を出す
// 廃棄処理: 元の例外を上書きしないため， stderr への記録にとどめる
void check_cuda(
    const CUresult result,
    const char* expression,
    const char* file,
    const int line);

[[nodiscard]]
bool report_cuda_cleanup_result(
    const CUresult result,
    const char* expression) noexcept;

} // namespace rainbow::detail

#define RAINBOW_CUDA_CHECK(expression) \
    do \
    { \
        ::rainbow::detail::check_cuda( \
            (expression), #expression, __FILE__, __LINE__ \
        ); \
    } while(false)
