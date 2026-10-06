# Include from tests/CMakeLists.txt. This test itself requires no GPU.
add_executable(rainbow_patch_failure_report_cpu_tests
    test_patch_failure_report_cpu.cpp
    "${PROJECT_SOURCE_DIR}/src/patch_failure_report.cpp"
)
target_include_directories(rainbow_patch_failure_report_cpu_tests PRIVATE
    "${PROJECT_SOURCE_DIR}/include")
if(MSVC)
    target_compile_options(rainbow_patch_failure_report_cpu_tests PRIVATE
        /utf-8 /W4 /Zc:__cplusplus /fp:strict)
    target_compile_definitions(rainbow_patch_failure_report_cpu_tests PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(rainbow_patch_failure_report_cpu_tests PRIVATE -fno-fast-math -ffp-contract=off)
endif()
add_test(NAME rainbow_patch_failure_report_cpu COMMAND rainbow_patch_failure_report_cpu_tests)
set_tests_properties(rainbow_patch_failure_report_cpu PROPERTIES LABELS "cpu;diagnostics;numerics")

add_executable(rainbow_patch_failure_capture_gpu_tests test_patch_failure_capture_gpu.cpp)
target_link_libraries(rainbow_patch_failure_capture_gpu_tests PRIVATE rainbow_drop_trace_host)
add_dependencies(rainbow_patch_failure_capture_gpu_tests
    rainbow_raindrop_trace_module rainbow_patch_build_module
    rainbow_patch_query_module rainbow_patch_optics_module)
if(MSVC)
    target_compile_options(rainbow_patch_failure_capture_gpu_tests PRIVATE /utf-8 /W4 /Zc:__cplusplus /fp:precise)
endif()
add_test(NAME rainbow_patch_failure_capture_gpu COMMAND rainbow_patch_failure_capture_gpu_tests
    "$<TARGET_OBJECTS:rainbow_raindrop_trace_module>"
    "$<TARGET_OBJECTS:rainbow_patch_build_module>"
    "$<TARGET_OBJECTS:rainbow_patch_query_module>"
    "$<TARGET_OBJECTS:rainbow_patch_optics_module>"
    COMMAND_EXPAND_LISTS)
set_tests_properties(rainbow_patch_failure_capture_gpu PROPERTIES LABELS "gpu;cuda;optix;diagnostics;integration")
