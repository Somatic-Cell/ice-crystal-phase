# Include from tests/CMakeLists.txt, AFTER patch_optics_tests.cmake.
# Reuse the production fatbin. No extra OptiX pipeline/module or CUDA runtime.
add_executable(rainbow_wave_optics_cpu_tests test_wave_optics_cpu.cpp)
target_include_directories(rainbow_wave_optics_cpu_tests PRIVATE "${PROJECT_SOURCE_DIR}/include")
add_executable(rainbow_wave_optics_cuda_tests test_wave_optics_cuda.cpp)
target_link_libraries(rainbow_wave_optics_cuda_tests PRIVATE rainbow_drop_trace_host)
add_dependencies(rainbow_wave_optics_cuda_tests
    rainbow_patch_optics_module rainbow_patch_build_module rainbow_patch_query_module)
foreach(wave_test IN ITEMS rainbow_wave_optics_cpu_tests rainbow_wave_optics_cuda_tests)
    if(MSVC)
        target_compile_options(${wave_test} PRIVATE /utf-8 /W4 /Zc:__cplusplus /fp:precise)
        target_compile_definitions(${wave_test} PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${wave_test} PRIVATE -fno-fast-math -ffp-contract=off)
    endif()
endforeach()
add_test(NAME rainbow_wave_optics_cpu COMMAND rainbow_wave_optics_cpu_tests)
add_test(NAME rainbow_wave_optics_cuda COMMAND rainbow_wave_optics_cuda_tests
    "$<TARGET_OBJECTS:rainbow_patch_optics_module>"
    "$<TARGET_OBJECTS:rainbow_patch_build_module>"
    "$<TARGET_OBJECTS:rainbow_patch_query_module>"
    COMMAND_EXPAND_LISTS)
set_tests_properties(rainbow_wave_optics_cpu PROPERTIES LABELS "cpu;optics;focal;diffraction")
set_tests_properties(rainbow_wave_optics_cuda PROPERTIES LABELS "gpu;cuda;optix;optics;integration")
