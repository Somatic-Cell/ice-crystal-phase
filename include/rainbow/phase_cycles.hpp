#pragma once

#include <rainbow/complex32.hpp>
#include <rainbow/host_device.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace rainbow
{

struct PhaseCycles
{
    std::int32_t turns = 0;
    float fraction = 0.0f;

    [[nodiscard]]
    HOST_DEVICE constexpr bool is_valid() const noexcept
    {
        // NaN だった場合，比較によって false が返る
        return fraction >= 0.0f && fraction < 1.0f;
    }

    // 整数部分と残りの周期数を別々に渡す
    [[nodiscard]]
    HOST_DEVICE bool try_set_parts(
        const std::int32_t whole_turns,
        const float residual_cycles
    ) noexcept
    {
        if(!(residual_cycles > -0x1p32f && residual_cycles < 0x1p32f)){
            return false;
        }
        return normalize_and_store(
            static_cast<std::int64_t>(whole_turns), residual_cycles, 0.0f
        );
    }

    [[nodiscard]]
    HOST_DEVICE bool try_add(const PhaseCycles rhs) noexcept
    {
        if(!is_valid() || !rhs.is_valid())
        {
            return false;
        }
        float sum = 0.0f;
        float error = 0.0f;
        two_sum(fraction, rhs.fraction, sum, error);
        return normalize_and_store(
            static_cast<std::int64_t>(turns) + static_cast<std::int64_t>(rhs.turns),
            sum,
            error
        );
    }

    [[nodiscard]]
    HOST_DEVICE bool try_subtract(const PhaseCycles rhs) noexcept
    {
        if(!is_valid() || !rhs.is_valid())
        {
            return false;
        }
        float difference = 0.0f;
        float error = 0.0f;
        two_sum(fraction, -rhs.fraction, difference, error);
        return normalize_and_store(
            static_cast<std::int64_t>(turns) - static_cast<std::int64_t>(rhs.turns),
            difference, error);
    }

    // exp(+i * 2*pi*q) を FP32 で評価する．整数周期を使わないのはこの最終段階だけ．
    // 光路・焦線の位相をいつ電場に適用するかは呼び出し側の責務である．
    // invalid な PhaseCycles の場合は false を返し，output を変更しない．
    [[nodiscard]]
    HOST_DEVICE bool try_unit_phasor(Complex32& output) const noexcept
    {
        if(!is_valid())
        {
            return false;
        }

        // 正確に表現できる四分の一周期は，三角関数の丸めを介さずに返す．
        if(fraction == 0.0f)  { output = { 1.0f,  0.0f}; return true; }
        if(fraction == 0.25f) { output = { 0.0f,  1.0f}; return true; }
        if(fraction == 0.5f)  { output = {-1.0f,  0.0f}; return true; }
        if(fraction == 0.75f) { output = { 0.0f, -1.0f}; return true; }

        // 大きな turns + fraction に 2*pi を掛けない．三角関数への引数は [-pi, pi]．
        // fraction >= 0.5 の 1 の減算は，この区間では正確に表現できる．
        const float reduced_cycles = fraction > 0.5f ? fraction - 1.0f : fraction;
        constexpr float two_pi = 0x1.921fb6p+2f;
        const float angle_rad = two_pi * reduced_cycles;
        output = {::cosf(angle_rad), ::sinf(angle_rad)};
        return true;
    }

private:
    HOST_DEVICE
    static void two_sum(
        const float lhs, const float rhs, float&sum, float& error
    ) noexcept
    {
        sum = lhs + rhs;
        const float rhs_virtual = sum - lhs;
        const float lhs_virtual = sum - rhs_virtual;
        error = (lhs - lhs_virtual) + (rhs - rhs_virtual);
    }

    [[nodiscard]]
    HOST_DEVICE bool normalize_and_store(
        std::int64_t whole_turns, const float hi, const float lo) noexcept
    {
        const float integral = ::floorf(hi);
        whole_turns += static_cast<std::int64_t>(integral);
        float new_fraction = (hi - integral) + lo;

        if(new_fraction < 0.0f)
        {
            --whole_turns;
            new_fraction += 1.0f;
        }
        // 負の微小値に 1 を足した結果が，float の丸めで 1 になる場合もここで扱う．
        if(new_fraction >= 1.0f)
        {
            ++whole_turns;
            new_fraction -= 1.0f;
        }

        if(whole_turns < -2147483648LL || whole_turns > 2147483647LL)
        {
            return false;
        }

        // 状態の更新は成功時だけ．小数部の -0 は +0 に正規化する．
        turns = static_cast<std::int32_t>(whole_turns);
        fraction = new_fraction == 0.0f ? 0.0f : new_fraction;
        return true;
    }


};

static_assert(std::numeric_limits<float>::is_iec559);
static_assert(std::numeric_limits<float>::digits == 24);
static_assert(sizeof(std::int32_t) == 4);
static_assert(std::is_standard_layout_v<PhaseCycles>);
static_assert(std::is_trivially_copyable_v<PhaseCycles>);
static_assert(sizeof(PhaseCycles) == 8);
static_assert(alignof(PhaseCycles) == 4);
static_assert(offsetof(PhaseCycles, turns) == 0);
static_assert(offsetof(PhaseCycles, fraction) == 4);

}