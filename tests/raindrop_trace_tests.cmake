# tests/CMakeLists.txt の末尾から include する．既存の 5 テストを削除・改名しない．
add_executable(rainbow_raindrop_cpu_tests test_raindrop_cpu.cpp)
target_include_directories(rainbow_raindrop_cpu_tests PRIVATE "${PROJECT_SOURCE_DIR}/include")
add_executable(rainbow_raindrop_optix_tests test_raindrop_optix.cpp)
target_link_libraries(rainbow_raindrop_optix_tests PRIVATE rainbow_drop_trace_host)
add_dependencies(rainbow_raindrop_optix_tests rainbow_raindrop_trace_module)
if(MSVC)
    foreach(target IN ITEMS rainbow_raindrop_cpu_tests rainbow_raindrop_optix_tests)
        target_compile_options(${target} PRIVATE /utf-8 /W4 /Zc:__cplusplus /fp:precise)
        target_compile_definitions(${target} PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
    endforeach()
else()
    foreach(target IN ITEMS rainbow_raindrop_cpu_tests rainbow_raindrop_optix_tests)
        target_compile_options(${target} PRIVATE -ffp-contract=off -fno-fast-math)
    endforeach()
endif()
add_test(NAME rainbow_raindrop_cpu COMMAND rainbow_raindrop_cpu_tests)
add_test(NAME rainbow_raindrop_optix
    COMMAND rainbow_raindrop_optix_tests "$<TARGET_OBJECTS:rainbow_raindrop_trace_module>"
    COMMAND_EXPAND_LISTS
)
set_tests_properties(rainbow_raindrop_cpu PROPERTIES LABELS "cpu;raindrop;physics")
set_tests_properties(rainbow_raindrop_optix PROPERTIES LABELS "gpu;optix;raindrop;physics")
