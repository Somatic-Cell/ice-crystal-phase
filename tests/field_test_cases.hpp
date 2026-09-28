#pragma once

// ホスト専用の入力生成・独立参照・検証．CUDA kernel からは include しない．
#include "field_test_data.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace rainbow::tests
{

[[nodiscard]]
inline std::vector<FieldTestInput> make_field_test_inputs()
{
    const Field32 horizontal{{1.0f, 0.0f}, {0.0f, 0.0f}};
    const Field32 vertical{{0.0f, 0.0f}, {1.0f, 0.0f}};
    const Field32 negative_horizontal{{-1.0f, 0.0f}, {0.0f, 0.0f}};
    std::vector<FieldTestInput> inputs;
    inputs.reserve(1027);

    // 先頭3例は，合成電場の二乗ノルムが解析的に 4 / 0 / 2 になる．
    inputs.push_back({horizontal, horizontal, {0.0f, 1.0f}, 2.0f, 0.0f, 1.0f});
    inputs.push_back({horizontal, negative_horizontal});
    inputs.push_back({horizontal, vertical});
    inputs.push_back({});
    inputs.push_back({{{1.0f, 2.0f}, {-3.0f, 0.5f}},
                      {{-2.0f, 1.0f}, {1.5f, -0.25f}},
                      {0.0f, -1.0f}, -0.5f, -1.0f, 0.0f});
    inputs.push_back({{{1.0e-4f, -2.0e-4f}, {3.0e-4f, 1.0e-4f}},
                      {{-1.0e-4f, 2.0e-4f}, {0.0f, 0.0f}},
                      {0.6f, 0.8f}, 0.125f, 0.6f, 0.8f});
    inputs.push_back({{{1.0000001f, 1.0f}, {1.0f, -1.0000001f}},
                      {{1.0f, 1.0000001f}, {-1.0f, 1.0f}},
                      {0.6f, 0.8f}, 0.0f, 0.6f, 0.8f});

    // 固定 seed + 整数からの変換により，標準ライブラリ固有の実数分布を避ける．
    std::mt19937 random(0x6ac19435u);
    const auto next_scalar = [&random]() -> float
    {
        const int centered = static_cast<int>(random() & 0xffffu) - 32768;
        return static_cast<float>(centered) / 8192.0f; // [-4,4) の有限値．
    };
    while(inputs.size() < 1027)
    {
        FieldTestInput input;
        input.left = {{next_scalar(), next_scalar()}, {next_scalar(), next_scalar()}};
        input.right = {{next_scalar(), next_scalar()}, {next_scalar(), next_scalar()}};
        input.factor = {0.5f * next_scalar(), 0.5f * next_scalar()};
        input.scale = 0.5f * next_scalar();
        const double angle_rad = static_cast<double>(next_scalar());
        input.cosine = static_cast<float>(std::cos(angle_rad));
        input.sine = static_cast<float>(std::sin(angle_rad));
        inputs.push_back(input);
    }
    // 1027 は 256 の倍数ではない．CUDA の末尾境界チェックも検査する．
    return inputs;
}

class FieldTestVerifier final
{
public:
    void verify(
        const std::span<const FieldTestInput> inputs,
        const std::span<const FieldTestResult> results)
    {
        if(inputs.size() != results.size() || inputs.size() < 3)
        {
            throw std::runtime_error("Field test input/output sizes are invalid.");
        }
        for(std::size_t index = 0; index < inputs.size(); ++index)
        {
            verify_one(inputs[index], results[index], index);
        }
        // 単純な解析例は float でも正確に表現できるため，許容誤差ではなく一致を確認．
        if(results[0].superposed_squared_norm != 4.0f
           || results[1].superposed_squared_norm != 0.0f
           || results[2].superposed_squared_norm != 2.0f)
        {
            throw std::runtime_error("Analytic coherent-sum checks (4, 0, 2) failed.");
        }
    }

    void print_summary(const std::string_view backend, const std::size_t case_count) const
    {
        std::cout << "Field arithmetic " << backend << ": passed ("
                  << case_count << " cases, " << scalar_check_count_ << " scalar checks)\n"
                  << "  maximum absolute error: " << std::setprecision(9)
                  << maximum_absolute_error_ << '\n';
    }

private:
    using ReferenceComplex = std::complex<double>;

    [[nodiscard]]
    static ReferenceComplex reference(const Complex32 value)
    {
        // 入力 float 自体の量子化は誤差比較に含めず，その値を FP64 に正確に変換．
        return {static_cast<double>(value.real), static_cast<double>(value.imag)};
    }

    void check_scalar(
        const float actual,
        const double expected,
        const double expression_scale,
        const std::size_t case_index,
        const std::string_view operation)
    {
        // 結果が0に近い減算にも対応するため，結果だけでなく入力項の大きさも使う．
        // 16 epsilon は，ここで扱う短い演算列と下記の有限入力集合の検証基準．
        // 任意の入力・任意の長い coherent sum の誤差保証ではない．
        const double scale = (std::max)(std::abs(expected), expression_scale);
        const double tolerance =
            16.0 * static_cast<double>(std::numeric_limits<float>::epsilon()) * scale
            + 1.0e-30;
        const double error = std::abs(static_cast<double>(actual) - expected);
        if(!std::isfinite(actual) || !std::isfinite(expected) || error > tolerance)
        {
            std::ostringstream message;
            message << std::setprecision(17) << "Field test case " << case_index
                    << ", " << operation << ": got " << actual
                    << ", expected " << expected << ", tolerance " << tolerance;
            throw std::runtime_error(message.str());
        }
        maximum_absolute_error_ = (std::max)(maximum_absolute_error_, error);
        ++scalar_check_count_;
    }

    void check_complex(
        const Complex32 actual,
        const ReferenceComplex expected,
        const double scale,
        const std::size_t index,
        const std::string_view operation)
    {
        check_scalar(actual.real, expected.real(), scale, index, operation);
        check_scalar(actual.imag, expected.imag(), scale, index, operation);
    }

    void verify_one(
        const FieldTestInput& input,
        const FieldTestResult& result,
        const std::size_t index)
    {
        if(result.completed_tag != (0x63af42d1u ^ static_cast<std::uint32_t>(index))
           || result.complex_byte_size != 8 || result.field_byte_size != 16)
        {
            throw std::runtime_error("Field test output tag / device ABI mismatch.");
        }

        // 参照演算は Complex32 / Field32 のメンバを一切呼ばない．
        const ReferenceComplex ax = reference(input.left.x);
        const ReferenceComplex ay = reference(input.left.y);
        const ReferenceComplex bx = reference(input.right.x);
        const ReferenceComplex by = reference(input.right.y);
        const ReferenceComplex factor = reference(input.factor);
        const double scale = static_cast<double>(input.scale);
        const double c = static_cast<double>(input.cosine);
        const double s = static_cast<double>(input.sine);
        const double ax_abs = std::abs(ax);
        const double ay_abs = std::abs(ay);
        const double bx_abs = std::abs(bx);
        const double by_abs = std::abs(by);

        check_complex(result.complex_sum, ax + bx, ax_abs + bx_abs, index, "complex sum");
        check_complex(result.complex_difference, ax - bx, ax_abs + bx_abs,
                      index, "complex difference");
        check_complex(result.complex_product, ax * bx, ax_abs * bx_abs,
                      index, "complex product");
        check_complex(result.complex_conjugate, std::conj(ax), ax_abs,
                      index, "complex conjugate");
        check_complex(result.complex_scaled, ax * scale, ax_abs * std::abs(scale),
                      index, "complex scale");
        check_complex(result.complex_square, ax * ax, ax_abs * ax_abs,
                      index, "complex self-multiply");
        check_scalar(result.complex_squared_norm, std::norm(ax), std::norm(ax),
                     index, "complex squared norm");
        check_complex(result.field_sum.x, ax + bx, ax_abs + bx_abs, index, "field sum.x");
        check_complex(result.field_sum.y, ay + by, ay_abs + by_abs, index, "field sum.y");
        check_complex(result.field_difference.x, ax - bx, ax_abs + bx_abs,
                      index, "field difference.x");
        check_complex(result.field_difference.y, ay - by, ay_abs + by_abs,
                      index, "field difference.y");
        check_complex(result.field_factored.x, ax * factor, ax_abs * std::abs(factor),
                      index, "field factor.x");
        check_complex(result.field_factored.y, ay * factor, ay_abs * std::abs(factor),
                      index, "field factor.y");
        check_complex(result.field_scaled.x, ax * scale, ax_abs * std::abs(scale),
                      index, "field scale.x");
        check_complex(result.field_scaled.y, ay * scale, ay_abs * std::abs(scale),
                      index, "field scale.y");
        check_complex(result.field_rotated.x, c * ax + s * ay,
                      std::abs(c) * ax_abs + std::abs(s) * ay_abs, index, "basis rotation.x");
        check_complex(result.field_rotated.y, -s * ax + c * ay,
                      std::abs(s) * ax_abs + std::abs(c) * ay_abs, index, "basis rotation.y");
        const double left_norm = std::norm(ax) + std::norm(ay);
        const double right_norm = std::norm(bx) + std::norm(by);
        check_scalar(result.field_squared_norm, left_norm, left_norm,
                     index, "field squared norm");
        check_scalar(result.superposed_squared_norm, std::norm(ax + bx) + std::norm(ay + by),
                     4.0 * (left_norm + right_norm), index, "coherent sum squared norm");

        // 結果の成分を double 化して，基底回転によるノルムの保存も独立に確認する．
        // c,s の float 化による単位円からの微小な誤差も，この検査の許容幅に含める．
        const double rotated_norm = std::norm(reference(result.field_rotated.x))
                                  + std::norm(reference(result.field_rotated.y));
        check_scalar(static_cast<float>(rotated_norm), left_norm,
                     2.0 * left_norm, index, "rotation preserves squared norm");
    }

    std::size_t scalar_check_count_ = 0;
    double maximum_absolute_error_ = 0.0;
};

} // namespace rainbow::tests
