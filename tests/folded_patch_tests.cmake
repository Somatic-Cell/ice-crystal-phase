# Include from tests/CMakeLists.txt after existing optical tests.
# Production fatbin and host library are reused; no new CUDA Runtime API.
add_executable(rainbow_folded_patch_cpu_tests test_folded_patches_cpu.cpp)
target_include_directories(rainbow_folded_patch_cpu_tests PRIVATE "${PROJECT_SOURCE_DIR}/include")
add_executable(rainbow_folded_patch_cuda_tests test_folded_patches_cuda.cpp)
target_link_libraries(rainbow_folded_patch_cuda_tests PRIVATE rainbow_drop_trace_host)
add_dependencies(rainbow_folded_patch_cuda_tests
    rainbow_patch_optics_module rainbow_patch_build_module rainbow_patch_query_module)
foreach(fold_target IN ITEMS rainbow_folded_patch_cpu_tests rainbow_folded_patch_cuda_tests)
    target_compile_features(${fold_target} PRIVATE cxx_std_20)
    if(MSVC)
        target_compile_options(${fold_target} PRIVATE /utf-8 /W4 /Zc:__cplusplus /fp:precise)
        target_compile_definitions(${fold_target} PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${fold_target} PRIVATE -fno-fast-math -ffp-contract=off)
    endif()
endforeach()
add_test(NAME rainbow_folded_patch_cpu COMMAND rainbow_folded_patch_cpu_tests)
add_test(NAME rainbow_folded_patch_cuda COMMAND rainbow_folded_patch_cuda_tests
    "$<TARGET_OBJECTS:rainbow_patch_optics_module>"
    "$<TARGET_OBJECTS:rainbow_patch_build_module>"
    "$<TARGET_OBJECTS:rainbow_patch_query_module>"
    COMMAND_EXPAND_LISTS)
set_tests_properties(rainbow_folded_patch_cpu PROPERTIES LABELS "cpu;optics;folded;regression")
set_tests_properties(rainbow_folded_patch_cuda PROPERTIES LABELS "gpu;cuda;optix;optics;integration")
