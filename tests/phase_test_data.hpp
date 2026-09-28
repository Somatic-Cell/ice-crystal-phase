#pragma once

#include <rainbow/phase_cycles.hpp>

#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace rainbow::tests
{

enum class PhaseTestOperation : std::uint32_t
{
    SetParts,
    Add,
    Subtract,
    UnitPhasor,
    RepeatedAdd
};

struct PhaseTestInput
{
    PhaseCycles initial = {};
    PhaseCycles operand = {};
    PhaseTestOperation operation = PhaseTestOperation::Add;
    std::uint32_t repetition_count = 1;
};

struct PhaseTestResult
{
    PhaseCycles phase = {};
    Complex32 phasor = {17.0f, -23.0f}; // 未変更の確認用 sentinel
    std::uint32_t succeeded = 0;
    std::uint32_t completed_tag = 0;
};

// 本番メンバを CPU/CUDA から同じ形で呼ぶだけ．独立参照計算は別ファイルにある．
[[nodiscard]]
HOST_DEVICE inline PhaseTestResult evaluate_phase_test(
    const PhaseTestInput& input, const std::uint32_t case_index) noexcept
{
    PhaseTestResult result = {};
    result.phase = input.initial;
    bool succeeded = false;
    switch(input.operation)
    {
    case PhaseTestOperation::SetParts:
        succeeded = result.phase.try_set_parts(input.operand.turns, input.operand.fraction);
        break;
    case PhaseTestOperation::Add:
        succeeded = result.phase.try_add(input.operand);
        break;
    case PhaseTestOperation::Subtract:
        succeeded = result.phase.try_subtract(input.operand);
        break;
    case PhaseTestOperation::UnitPhasor:
        succeeded = result.phase.try_unit_phasor(result.phasor);
        break;
    case PhaseTestOperation::RepeatedAdd:
        succeeded = true;
        for(std::uint32_t index = 0; index < input.repetition_count && succeeded; ++index)
        {
            succeeded = result.phase.try_add(input.operand);
        }
        break;
    }
    result.succeeded = succeeded ? 1u : 0u;
    result.completed_tag = 0x38d95a71u ^ case_index;
    return result;
}

static_assert(std::is_standard_layout_v<PhaseTestInput>);
static_assert(std::is_trivially_copyable_v<PhaseTestInput>);
static_assert(sizeof(PhaseTestInput) == 24);
static_assert(offsetof(PhaseTestInput, operation) == 16);
static_assert(std::is_standard_layout_v<PhaseTestResult>);
static_assert(std::is_trivially_copyable_v<PhaseTestResult>);
static_assert(sizeof(PhaseTestResult) == 24);
static_assert(offsetof(PhaseTestResult, completed_tag) == 20);

} // namespace rainbow::tests
