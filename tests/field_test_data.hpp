#pragma once

#include <rainbow/field32.hpp>

#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace rainbow::tests
{

// テスト用の転送 ABI．本番の launch parameter / OutgoingVertex ではない．
struct FieldTestInput
{
    Field32 left = {};
    Field32 right = {};
    Complex32 factor = {1.0f, 0.0f};
    float scale = 1.0f;
    float cosine = 1.0f;
    float sine = 0.0f;
};

struct FieldTestResult
{
    Complex32 complex_sum = {};
    Complex32 complex_difference = {};
    Complex32 complex_product = {};
    Complex32 complex_conjugate = {};
    Complex32 complex_scaled = {};
    Complex32 complex_square = {};
    float complex_squared_norm = 0.0f;
    Field32 field_sum = {};
    Field32 field_difference = {};
    Field32 field_factored = {};
    Field32 field_scaled = {};
    Field32 field_rotated = {};
    float field_squared_norm = 0.0f;
    float superposed_squared_norm = 0.0f;
    std::uint32_t complex_byte_size = 0;
    std::uint32_t field_byte_size = 0;
    std::uint32_t completed_tag = 0;
};

// CPU と CUDA が同じ「本番演算」を呼ぶための薄いテストアダプタ．
// 正解の計算ではない．参照値は field_test_cases.hpp 内で独立に求める．
[[nodiscard]]
HOST_DEVICE inline FieldTestResult evaluate_field_test(
    const FieldTestInput& input,
    const std::uint32_t case_index) noexcept
{
    FieldTestResult result = {};
    result.complex_sum = input.left.x + input.right.x;
    result.complex_difference = input.left.x - input.right.x;
    result.complex_product = input.left.x * input.right.x;
    result.complex_conjugate = input.left.x.conjugate();
    result.complex_scaled = input.left.x * input.scale;

    Complex32 square = input.left.x;
    square *= square; // 自己代入的な演算の検査．
    result.complex_square = square;
    result.complex_squared_norm = input.left.x.squared_norm();

    result.field_sum = input.left;
    result.field_sum += input.right;
    result.field_difference = input.left;
    result.field_difference -= input.right;
    result.field_factored = input.left;
    result.field_factored *= input.factor;
    result.field_scaled = input.left * input.scale;
    result.field_rotated = input.left.rotated_basis(input.cosine, input.sine);
    result.field_squared_norm = input.left.squared_norm();
    result.superposed_squared_norm = result.field_sum.squared_norm();

    // device 側からもレイアウトと全ケースの書き込みを確認する．
    result.complex_byte_size = static_cast<std::uint32_t>(sizeof(Complex32));
    result.field_byte_size = static_cast<std::uint32_t>(sizeof(Field32));
    result.completed_tag = 0x63af42d1u ^ case_index;
    return result;
}

static_assert(std::is_standard_layout_v<FieldTestInput>);
static_assert(std::is_trivially_copyable_v<FieldTestInput>);
static_assert(sizeof(FieldTestInput) == 52);
static_assert(offsetof(FieldTestInput, factor) == 32);
static_assert(offsetof(FieldTestInput, sine) == 48);
static_assert(std::is_standard_layout_v<FieldTestResult>);
static_assert(std::is_trivially_copyable_v<FieldTestResult>);
static_assert(sizeof(FieldTestResult) == 152);
static_assert(offsetof(FieldTestResult, field_sum) == 52);
static_assert(offsetof(FieldTestResult, completed_tag) == 148);

} // namespace rainbow::tests
