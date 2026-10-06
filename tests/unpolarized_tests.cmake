# Include once at the end of tests/CMakeLists.txt. Existing tests remain enabled.
add_executable(rainbow_unpolarized_cpu_tests test_unpolarized_cpu.cpp
    "${PROJECT_SOURCE_DIR}/src/patch_failure_report.cpp")
target_include_directories(rainbow_unpolarized_cpu_tests PRIVATE "${PROJECT_SOURCE_DIR}/include")
add_executable(rainbow_unpolarized_cuda_tests test_unpolarized_cuda.cpp)
target_link_libraries(rainbow_unpolarized_cuda_tests PRIVATE rainbow_drop_trace_host)
add_dependencies(rainbow_unpolarized_cuda_tests rainbow_raindrop_trace_module
    rainbow_patch_build_module rainbow_patch_query_module rainbow_patch_optics_module)
foreach(t IN ITEMS rainbow_unpolarized_cpu_tests rainbow_unpolarized_cuda_tests)
    target_compile_features(${t} PRIVATE cxx_std_20)
    if(MSVC)
        target_compile_options(${t} PRIVATE /utf-8 /W4 /Zc:__cplusplus /fp:strict)
        target_compile_definitions(${t} PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${t} PRIVATE -fno-fast-math -ffp-contract=off)
    endif()
endforeach()
add_test(NAME rainbow_unpolarized_cpu COMMAND rainbow_unpolarized_cpu_tests)
add_test(NAME rainbow_unpolarized_cuda COMMAND rainbow_unpolarized_cuda_tests
    "$<TARGET_OBJECTS:rainbow_raindrop_trace_module>"
    "$<TARGET_OBJECTS:rainbow_patch_build_module>"
    "$<TARGET_OBJECTS:rainbow_patch_query_module>"
    "$<TARGET_OBJECTS:rainbow_patch_optics_module>"
    COMMAND_EXPAND_LISTS)
set_tests_properties(rainbow_unpolarized_cpu PROPERTIES LABELS "cpu;optics;unpolarized;regression")
set_tests_properties(rainbow_unpolarized_cuda PROPERTIES LABELS "gpu;cuda;optix;unpolarized;integration")
