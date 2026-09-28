#include "field_test_data.hpp"

#include <cstdint>

// 通常 CUDA の数値テスト用 kernel．雨粒追跡・パッチ構築は行わない．
// kernel は入力と出力を対応付けるだけで，演算本体を複製しない．
extern "C" __global__
void evaluate_field_arithmetic(
    const rainbow::tests::FieldTestInput* inputs,
    rainbow::tests::FieldTestResult* outputs,
    const std::uint32_t count)
{
    const std::uint32_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index < count)
    {
        outputs[index] = rainbow::tests::evaluate_field_test(inputs[index], index);
    }
}
