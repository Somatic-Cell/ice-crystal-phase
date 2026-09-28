#pragma once

#include "phase_test_data.hpp"

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>
#include <span>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace rainbow::tests
{

[[nodiscard]]
inline std::vector<PhaseTestInput> make_phase_test_inputs()
{
    using Operation = PhaseTestOperation;
    constexpr std::int32_t lower = (std::numeric_limits<std::int32_t>::min)();
    constexpr std::int32_t upper = (std::numeric_limits<std::int32_t>::max)();
    const float infinity = std::numeric_limits<float>::infinity();
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float below_one = std::nextafter(1.0f, 0.0f);
    const float above_quarter = std::nextafter(0.25f, infinity);
    const float below_quarter = std::nextafter(0.25f, 0.0f);
    std::vector<PhaseTestInput> inputs;

    const auto add_case = [&](const PhaseCycles initial, const PhaseCycles operand,
                              const Operation operation, const std::uint32_t count = 1)
    {
        inputs.push_back({initial, operand, operation, count});
    };

    // 非正規化の小数部，符号付きゼロ，保存範囲端，非有限値．
    for(const float residual : {0.0f, -0.0f, 0.25f, 1.0f, 1.25f, -0.25f, -1.0f,
                               -1.25f, below_one, 0x1p-75f, -0x1p-75f,
                               infinity, -infinity, nan})
    {
        add_case({37, 0.375f}, {10, residual}, Operation::SetParts);
    }
    add_case({37, 0.375f}, {lower, 0x1p31f}, Operation::SetParts);
    add_case({37, 0.375f}, {upper, -0x1p31f}, Operation::SetParts);
    add_case({37, 0.375f}, {0, 0x1p32f}, Operation::SetParts);
    add_case({37, 0.375f}, {0, -0x1p32f}, Operation::SetParts);
    add_case({37, 0.375f}, {upper, 1.0f}, Operation::SetParts);
    add_case({37, 0.375f}, {lower, -0.25f}, Operation::SetParts);

    // 1.0f に丸められた和を単純に減算した実装では，小さな余りを失うケース．
    add_case({1000000000, 0.75f}, {0, above_quarter}, Operation::Add);
    add_case({1000000000, 0.75f}, {0, below_quarter}, Operation::Add);
    add_case({1000000000, 0.75f}, {0, 0.25f}, Operation::Add);
    add_case({1000000000, 0.25f}, {0, 0.75f}, Operation::Subtract);
    add_case({1000000001, 0.02f}, {1000000000, 0.98f}, Operation::Subtract);
    add_case({upper, 0.75f}, {0, 0.5f}, Operation::Add);
    add_case({lower, 0.25f}, {0, 0.5f}, Operation::Subtract);
    add_case({upper, 0.0f}, {lower, 0.0f}, Operation::Subtract);
    add_case({upper, 0.25f}, {lower, 0.5f}, Operation::Add);
    add_case({0, 0x1p-149f}, {0, 0x1p-149f}, Operation::Add);

    for(const float invalid : {-0.25f, 1.0f, infinity, nan})
    {
        add_case({7, 0.5f}, {0, invalid}, Operation::Add);
        add_case({7, invalid}, {0, 0.25f}, Operation::Subtract);
        add_case({7, invalid}, {}, Operation::UnitPhasor);
    }

    // 整数周期の大小・符号を変えても，位相因子は小数部だけで決まる．
    for(const std::int32_t turns : {0, 1000000000, lower, upper})
    {
        for(const float fraction : {0.0f, 0.25f, 0.5f, 0.75f, below_one, 0x1p-75f})
        {
            add_case({turns, fraction}, {}, Operation::UnitPhasor);
        }
    }
    // この二進分数の反復和は保存精度に収まる．巨大な全周期へ float で足す方式は失敗する．
    add_case({1000000000, 0.0f}, {0, 0x1p-20f}, Operation::RepeatedAdd, 65537);
    add_case({upper, 0.5f}, {0, 0.25f}, Operation::RepeatedAdd, 3);

    // 固定 seed．double 参照の加減算で二進分数を正確に扱える入力を生成する．
    std::mt19937 generator{0x92ac4071u};
    for(std::size_t index = 0; index < 512; ++index)
    {
        const auto random_turns = [&]() -> std::int32_t
        {
            return static_cast<std::int32_t>(generator() % 2000000001u) - 1000000000;
        };
        const auto random_fraction = [&]() -> float
        {
            return static_cast<float>(generator() & 0xffffffu) * 0x1p-24f;
        };
        const PhaseCycles lhs{random_turns(), random_fraction()};
        const PhaseCycles rhs{random_turns(), random_fraction()};
        add_case(lhs, rhs, Operation::Add);
        add_case(lhs, rhs, Operation::Subtract);
        add_case(lhs, {}, Operation::UnitPhasor);
    }
    // 256 の倍数で終了しないようにし，kernel の末尾の範囲チェックも検査する．
    if(inputs.size() % 256 == 0)
    {
        add_case({0, 0.5f}, {}, Operation::UnitPhasor);
    }
    return inputs;
}

class PhaseTestVerifier final
{
public:
    void verify(const std::span<const PhaseTestInput> inputs,
                const std::span<const PhaseTestResult> outputs)
    {
        if(inputs.size() != outputs.size())
        {
            throw std::runtime_error("Phase test input/output count mismatch.");
        }
        for(std::size_t index = 0; index < inputs.size(); ++index)
        {
            const auto& input = inputs[index];
            const auto& output = outputs[index];
            PhaseCycles expected = input.initial;
            bool success = false;
            switch(input.operation)
            {
            case PhaseTestOperation::SetParts:
                success = reference_set(input.operand.turns, input.operand.fraction, expected);
                break;
            case PhaseTestOperation::Add:
                success = reference_arithmetic(expected, input.operand, false);
                break;
            case PhaseTestOperation::Subtract:
                success = reference_arithmetic(expected, input.operand, true);
                break;
            case PhaseTestOperation::UnitPhasor:
                success = valid(expected);
                break;
            case PhaseTestOperation::RepeatedAdd:
                success = true;
                for(std::uint32_t step = 0; step < input.repetition_count && success; ++step)
                {
                    success = reference_arithmetic(expected, input.operand, false);
                }
                break;
            }

            require(output.completed_tag == (0x38d95a71u ^ static_cast<std::uint32_t>(index)),
                    index, "kernel completion tag");
            require(output.succeeded == (success ? 1u : 0u), index, "success/failure status");
            require(output.phase.turns == expected.turns, index, "integer turns");
            require(same_bits(output.phase.fraction, expected.fraction), index, "fraction bits");

            if(input.operation == PhaseTestOperation::UnitPhasor && success)
            {
                verify_phasor(input.initial.fraction, output.phasor, index);
            }
            else
            {
                require(same_bits(output.phasor.real, 17.0f)
                        && same_bits(output.phasor.imag, -23.0f), index, "unchanged phasor");
            }
        }
    }

    void print_summary(const char* backend, const std::size_t case_count) const
    {
        std::cout << "PhaseCycles " << backend << " test: passed (" << case_count
                  << " cases, " << comparison_count_ << " checks)\n";
    }

private:
    // 参照側では PhaseCycles の演算メンバも TwoSum も使わない．
    // 整数を int64，小数の和差を double で独立に計算し，保存時にだけ float に丸める．
    [[nodiscard]]
    static bool valid(const PhaseCycles value)
    {
        return std::isfinite(value.fraction)
            && value.fraction >= 0.0f && value.fraction < 1.0f;
    }

    [[nodiscard]]
    static bool reference_set(std::int64_t turns, const double residual, PhaseCycles& out)
    {
        if(!std::isfinite(residual) || residual <= -0x1p32 || residual >= 0x1p32)
        {
            return false;
        }
        const double integer_part = std::floor(residual);
        turns += static_cast<std::int64_t>(integer_part);
        float fraction = static_cast<float>(residual - integer_part);
        if(fraction == 1.0f)
        {
            ++turns;
            fraction = 0.0f;
        }
        if(turns < -2147483648LL || turns > 2147483647LL)
        {
            return false;
        }
        out = {static_cast<std::int32_t>(turns), fraction == 0.0f ? 0.0f : fraction};
        return true;
    }

    [[nodiscard]]
    static bool reference_arithmetic(PhaseCycles& lhs, const PhaseCycles rhs,
                                     const bool subtract)
    {
        if(!valid(lhs) || !valid(rhs))
        {
            return false;
        }
        const std::int64_t turns = subtract
            ? static_cast<std::int64_t>(lhs.turns) - rhs.turns
            : static_cast<std::int64_t>(lhs.turns) + rhs.turns;
        const double fraction = subtract
            ? static_cast<double>(lhs.fraction) - static_cast<double>(rhs.fraction)
            : static_cast<double>(lhs.fraction) + static_cast<double>(rhs.fraction);
        return reference_set(turns, fraction, lhs);
    }

    [[nodiscard]]
    static bool same_bits(const float lhs, const float rhs)
    {
        return std::bit_cast<std::uint32_t>(lhs) == std::bit_cast<std::uint32_t>(rhs);
    }

    void require(const bool condition, const std::size_t index, const char* label)
    {
        ++comparison_count_;
        if(!condition)
        {
            std::ostringstream message;
            message << "Phase test case " << index << " failed: " << label;
            throw std::runtime_error(message.str());
        }
    }

    void verify_phasor(const float fraction, const Complex32 value, const std::size_t index)
    {
        // 本番の float 2*pi と reduction を流用せず，double で参照を求める．
        const double angle = 2.0 * std::numbers::pi_v<double> * static_cast<double>(fraction);
        constexpr double component_tolerance = 4.0e-7;
        require(std::isfinite(value.real) && std::isfinite(value.imag), index, "finite phasor");
        require(std::abs(static_cast<double>(value.real) - std::cos(angle)) <= component_tolerance,
                index, "phasor real");
        require(std::abs(static_cast<double>(value.imag) - std::sin(angle)) <= component_tolerance,
                index, "phasor imag");
        // 任意の角度を正確に単位円上の float へ置ける，という保証ではない．
        const double norm = static_cast<double>(value.real) * value.real
                          + static_cast<double>(value.imag) * value.imag;
        require(std::abs(norm - 1.0) <= 5.0e-7, index, "unit phasor squared norm");
        if(fraction == 0.0f || fraction == 0.25f || fraction == 0.5f || fraction == 0.75f)
        {
            const Complex32 expected = fraction == 0.0f ? Complex32{1.0f, 0.0f}
                : fraction == 0.25f ? Complex32{0.0f, 1.0f}
                : fraction == 0.5f ? Complex32{-1.0f, 0.0f} : Complex32{0.0f, -1.0f};
            require(value.real == expected.real && value.imag == expected.imag,
                    index, "exact quarter-cycle phasor");
        }
    }

    std::size_t comparison_count_ = 0;
};

} // namespace rainbow::tests
