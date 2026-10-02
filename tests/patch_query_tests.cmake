# tests/CMakeLists.txt から読む．root から直接 include しない．
add_executable(rainbow_patch_query_cpu_tests test_patch_query_cpu.cpp)
target_include_directories(rainbow_patch_query_cpu_tests PRIVATE "${PROJECT_SOURCE_DIR}/include")
add_executable(rainbow_patch_query_optix_tests test_patch_query_optix.cpp)
target_link_libraries(rainbow_patch_query_optix_tests PRIVATE rainbow_drop_trace_host)
add_dependencies(rainbow_patch_query_optix_tests rainbow_patch_build_module rainbow_patch_query_module)
foreach(query_test IN ITEMS rainbow_patch_query_cpu_tests rainbow_patch_query_optix_tests)
    if(MSVC)
        target_compile_options(${query_test} PRIVATE /utf-8 /W4 /Zc:__cplusplus /fp:precise)
        target_compile_definitions(${query_test} PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
    else()
        target_compile_options(${query_test} PRIVATE -ffp-contract=off -fno-fast-math)
    endif()
endforeach()
add_test(NAME rainbow_patch_query_cpu COMMAND rainbow_patch_query_cpu_tests)
add_test(NAME rainbow_patch_query_optix
    COMMAND rainbow_patch_query_optix_tests
        "$<TARGET_OBJECTS:rainbow_patch_build_module>"
        "$<TARGET_OBJECTS:rainbow_patch_query_module>"
    COMMAND_EXPAND_LISTS
)
set_tests_properties(rainbow_patch_query_cpu PROPERTIES LABELS "cpu;patch;query;geometry")
set_tests_properties(rainbow_patch_query_optix PROPERTIES LABELS "gpu;optix;patch;query;geometry")
