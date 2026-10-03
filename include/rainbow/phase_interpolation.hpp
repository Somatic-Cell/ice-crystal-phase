#pragma once

#include <rainbow/phase_cycles.hpp>

#include <cmath>
#include <cstdint>
#include <limits>

namespace rainbow
{

// Bilinear interpolation of optical path measured in cycles.
//
// Important: do not interpolate exp(i 2 pi q) at the corners.  The paper
// interpolates the optical path first and applies propagation phase after the
// interpolation.  We retain q00 as an integer+fraction anchor and interpolate
// only path differences in double precision, avoiding the loss that would
// result from first forming a large float (turns + fraction).
struct PhaseInterpolation
{
    [[nodiscard]]
    HOST_DEVICE static bool try_bilinear(
        const PhaseCycles q00,
        const PhaseCycles q10,
        const PhaseCycles q01,
        const PhaseCycles q11,
        const float u,
        const float v,
        PhaseCycles& output) noexcept
    {
        if(!q00.is_valid() || !q10.is_valid() || !q01.is_valid() || !q11.is_valid()
           || !(u >= 0.0f && u <= 1.0f)
           || !(v >= 0.0f && v <= 1.0f))
        {
            return false;
        }

        // Preserve corners bit-for-bit.  Besides being useful numerically, this
        // makes the interpolation contract explicit at patch boundaries.
        if(u == 0.0f && v == 0.0f) { output = q00; return true; }
        if(u == 1.0f && v == 0.0f) { output = q10; return true; }
        if(u == 0.0f && v == 1.0f) { output = q01; return true; }
        if(u == 1.0f && v == 1.0f) { output = q11; return true; }

        // q(u,v) = q00 + u (q10-q00)
        //          + v [(q01-q00) + u{(q11-q01)-(q10-q00)}].
        //
        // The turn differences are formed in int64_t before conversion, so no
        // signed overflow occurs even when two valid int32_t endpoints are far
        // apart.  In the raindrop application neighbouring vertices differ by
        // only a small number of optical cycles; double therefore preserves far
        // more phase information than the final FP32 phasor evaluation needs.
        const double du0 = difference(q10, q00);
        const double du1 = difference(q11, q01);
        const double dv0 = difference(q01, q00);
        const double ud = static_cast<double>(u);
        const double vd = static_cast<double>(v);
        const double row_difference = ::fma(ud, du1 - du0, dv0);
        const double delta = ::fma(vd, row_difference, ud * du0);

        return add_delta(q00, delta, output);
    }

private:
    // Scalar constants, avoiding host-only numeric_limits calls in device code.
    static constexpr std::int64_t minimum_turns = -2147483647LL - 1LL;
    static constexpr std::int64_t maximum_turns = 2147483647LL;
    [[nodiscard]]
    HOST_DEVICE static double difference(
        const PhaseCycles lhs,
        const PhaseCycles rhs) noexcept
    {
        const std::int64_t whole =
            static_cast<std::int64_t>(lhs.turns)
            - static_cast<std::int64_t>(rhs.turns);
        return static_cast<double>(whole)
            + (static_cast<double>(lhs.fraction)
               - static_cast<double>(rhs.fraction));
    }

    [[nodiscard]]
    HOST_DEVICE static bool add_delta(
        const PhaseCycles anchor,
        const double delta,
        PhaseCycles& output) noexcept
    {
        // NaN fails both comparisons.  The bound is wider than every possible
        // difference of two int32_t turns, but still makes the conversion to
        // int64_t explicit and well-defined.
        constexpr double delta_limit = 0x1p33;
        if(!(delta > -delta_limit && delta < delta_limit))
        {
            return false;
        }

        const double integral_double = ::floor(delta);
        const auto integral = static_cast<std::int64_t>(integral_double);
        double fraction =
            static_cast<double>(anchor.fraction)
            + (delta - integral_double);
        std::int64_t turns =
            static_cast<std::int64_t>(anchor.turns) + integral;

        // delta-floor(delta) is in [0,1), while anchor.fraction is in [0,1).
        // One carry is therefore sufficient.  The lower branch only guards
        // against a possible one-ulp roundoff below zero.
        if(fraction < 0.0)
        {
            --turns;
            fraction += 1.0;
        }
        if(fraction >= 1.0)
        {
            ++turns;
            fraction -= 1.0;
        }

        if(turns < minimum_turns
           || turns > maximum_turns)
        {
            return false;
        }

        float stored_fraction = static_cast<float>(fraction);
        // A double value just below one can round to exactly 1 in float.
        if(stored_fraction >= 1.0f)
        {
            ++turns;
            stored_fraction = 0.0f;
            if(turns > maximum_turns)
            {
                return false;
            }
        }
        else if(stored_fraction < 0.0f)
        {
            --turns;
            stored_fraction += 1.0f;
            if(turns < minimum_turns)
            {
                return false;
            }
        }

        output.turns = static_cast<std::int32_t>(turns);
        output.fraction = stored_fraction == 0.0f ? 0.0f : stored_fraction;
        return output.is_valid();
    }
};

} // namespace rainbow
