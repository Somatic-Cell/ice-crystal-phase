#include "field_test_cases.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <span>
#include <vector>

int main()
{
    try
    {
        const auto inputs = rainbow::tests::make_field_test_inputs();
        std::vector<rainbow::tests::FieldTestResult> outputs(inputs.size());
        for(std::size_t index = 0; index < inputs.size(); ++index)
        {
            outputs[index] = rainbow::tests::evaluate_field_test(
                inputs[index], static_cast<std::uint32_t>(index));
        }
        rainbow::tests::FieldTestVerifier verifier;
        verifier.verify(
            std::span<const rainbow::tests::FieldTestInput>{inputs},
            std::span<const rainbow::tests::FieldTestResult>{outputs});
        verifier.print_summary("CPU", inputs.size());
        return EXIT_SUCCESS;
    }
    catch(const std::exception& exception)
    {
        std::cerr << "Error: " << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
