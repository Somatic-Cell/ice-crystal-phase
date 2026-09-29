# Include from tests/CMakeLists.txt. Do not include this directly from root.
add_executable(rainbow_patch_cpu_tests test_patch_cpu.cpp)
target_include_directories(rainbow_patch_cpu_tests PRIVATE "${PROJECT_SOURCE_DIR}/include")

add_executable(rainbow_patch_optix_tests test_patch_optix.cpp)
target_link_libraries(rainbow_patch_optix_tests PRIVATE rainbow_drop_trace_host)
add_dependencies(rainbow_patch_optix_tests rainbow_patch_build_module rainbow_raindrop_trace_module)

foreach(patch_test IN ITEMS rainbow_patch_cpu_tests rainbow_patch_optix_tests)
    if(MSVC)
        target_compile_options(${patch_test} PRIVATE /utf-8 /W4 /Zc:__cplusplus /fp:precise)
        target_compile_definitions(${patch_test} PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
    else()
        target_compile_options(${patch_test} PRIVATE -ffp-contract=off -fno-fast-math)
    endif()
endforeach()
add_test(NAME rainbow_patch_cpu COMMAND rainbow_patch_cpu_tests)
add_test(NAME rainbow_patch_optix
    COMMAND rainbow_patch_optix_tests
        "$<TARGET_OBJECTS:rainbow_patch_build_module>"
        "$<TARGET_OBJECTS:rainbow_raindrop_trace_module>"
    COMMAND_EXPAND_LISTS
)
set_tests_properties(rainbow_patch_cpu PROPERTIES LABELS "cpu;patch;geometry")
set_tests_properties(rainbow_patch_optix PROPERTIES LABELS "gpu;cuda;optix;patch;geometry")
