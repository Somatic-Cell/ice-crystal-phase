#include "field_test_cases.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <span>
#include <vector>
#include <type_traits>
#include <utility>

// 通常の複合代入と同様に，左辺自身への参照を返す契約を検査する．
// 値を返すと，(a += b) += c の二回目は a ではなく一時値を更新してしまう．
static_assert(std::is_same_v<
    decltype(std::declval<rainbow::Field32&>() += std::declval<rainbow::Field32>()),
    rainbow::Field32&>);
static_assert(std::is_same_v<
    decltype(std::declval<rainbow::Field32&>() -= std::declval<rainbow::Field32>()),
    rainbow::Field32&>);

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
