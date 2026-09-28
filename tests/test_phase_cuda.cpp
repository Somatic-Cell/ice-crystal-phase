#include "phase_test_cases.hpp"
#include "cuda_array_test.hpp"

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>

// Windows では CTest が渡す module パスをワイド文字で受け取る．
#if defined(_WIN32)
int wmain(const int argc, wchar_t* argv[])
#else
int main(const int argc, char* argv[])
#endif
{
    if(argc != 2)
    {
        std::cerr << "Usage: rainbow_phase_cuda_tests <phase_test.fatbin>\n";
        return EXIT_FAILURE;
    }
    try
    {
        rainbow::CudaContext cuda_context{0};
        rainbow::tests::CudaArrayTest<
            rainbow::tests::PhaseTestInput,
            rainbow::tests::PhaseTestResult> test(
                cuda_context, rainbow::tests::make_phase_test_inputs());
        test.run(std::filesystem::path{argv[1]}, "evaluate_phase_arithmetic");

        // 入力生成と数値検証はテスト固有．GPU 資源管理は field / phase で共有する．
        rainbow::tests::PhaseTestVerifier verifier;
        verifier.verify(test.inputs(), test.outputs());
        verifier.print_summary("CUDA", test.inputs().size());
        return EXIT_SUCCESS;
    }
    catch(const std::exception& exception)
    {
        std::cerr << "Error: " << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
