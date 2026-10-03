add_executable(
    rainbow_optical_interpolation_cpu_tests
    test_optical_interpolation_cpu.cpp
)

target_include_directories(
    rainbow_optical_interpolation_cpu_tests
    PRIVATE
        "${PROJECT_SOURCE_DIR}/include"
)

if(MSVC)
    target_compile_options(
        rainbow_optical_interpolation_cpu_tests
        PRIVATE
            /utf-8 /W4 /Zc:__cplusplus /fp:precise
    )
    target_compile_definitions(
        rainbow_optical_interpolation_cpu_tests
        PRIVATE
            NOMINMAX WIN32_LEAN_AND_MEAN
    )
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(
        rainbow_optical_interpolation_cpu_tests
        PRIVATE
            -fno-fast-math -ffp-contract=off
    )
endif()

add_test(
    NAME rainbow_optical_interpolation_cpu
    COMMAND rainbow_optical_interpolation_cpu_tests
)

set_tests_properties(
    rainbow_optical_interpolation_cpu
    PROPERTIES
        LABELS "cpu;numerics;optics;interpolation"
)
