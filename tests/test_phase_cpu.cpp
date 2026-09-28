#include "phase_test_cases.hpp"

#include <cfenv>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <span>
#include <vector>

int main()
{
    try
    {
        if(std::fegetround() != FE_TONEAREST)
        {
            throw std::runtime_error("Phase tests require round-to-nearest mode.");
        }
        const auto inputs = rainbow::tests::make_phase_test_inputs();
        std::vector<rainbow::tests::PhaseTestResult> outputs(inputs.size());
        for(std::size_t index = 0; index < inputs.size(); ++index)
        {
            outputs[index] = rainbow::tests::evaluate_phase_test(
                inputs[index], static_cast<std::uint32_t>(index));
        }
        rainbow::tests::PhaseTestVerifier verifier;
        verifier.verify(inputs, outputs);
        verifier.print_summary("CPU", inputs.size());
        return EXIT_SUCCESS;
    }
    catch(const std::exception& exception)
    {
        std::cerr << "Error: " << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
