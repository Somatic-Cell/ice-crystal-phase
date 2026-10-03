# Included from tests/CMakeLists.txt, not from the root directory.
add_executable(rainbow_patch_optics_cpu_tests test_patch_optics_cpu.cpp)
target_include_directories(rainbow_patch_optics_cpu_tests PRIVATE "${PROJECT_SOURCE_DIR}/include")
add_executable(rainbow_patch_optics_cuda_tests test_patch_optics_cuda.cpp)
target_link_libraries(rainbow_patch_optics_cuda_tests PRIVATE rainbow_drop_trace_host)
add_dependencies(rainbow_patch_optics_cuda_tests
    rainbow_patch_optics_module rainbow_patch_build_module rainbow_patch_query_module)
foreach(optical_test IN ITEMS rainbow_patch_optics_cpu_tests rainbow_patch_optics_cuda_tests)
    if(MSVC)
        target_compile_options(${optical_test} PRIVATE /utf-8 /W4 /Zc:__cplusplus /fp:precise)
        target_compile_definitions(${optical_test} PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${optical_test} PRIVATE -fno-fast-math -ffp-contract=off)
    endif()
endforeach()
add_test(NAME rainbow_patch_optics_cpu COMMAND rainbow_patch_optics_cpu_tests)
add_test(NAME rainbow_patch_optics_cuda COMMAND rainbow_patch_optics_cuda_tests
    "$<TARGET_OBJECTS:rainbow_patch_optics_module>"
    "$<TARGET_OBJECTS:rainbow_patch_build_module>"
    "$<TARGET_OBJECTS:rainbow_patch_query_module>"
    COMMAND_EXPAND_LISTS
)
set_tests_properties(rainbow_patch_optics_cpu PROPERTIES LABELS "cpu;optics;numerics")
set_tests_properties(rainbow_patch_optics_cuda PROPERTIES LABELS "gpu;cuda;optix;optics;integration")
