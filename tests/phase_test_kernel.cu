#include "phase_test_data.hpp"

#include <cstdint>

extern "C" __global__
void evaluate_phase_arithmetic(
    const rainbow::tests::PhaseTestInput* inputs,
    rainbow::tests::PhaseTestResult* outputs,
    const std::uint32_t count)
{
    const std::uint32_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index < count)
    {
        outputs[index] = rainbow::tests::evaluate_phase_test(inputs[index], index);
    }
}
