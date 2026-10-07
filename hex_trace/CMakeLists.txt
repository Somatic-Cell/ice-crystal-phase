cmake_minimum_required(VERSION 3.27)

# Also usable as a standalone build: -S hex_trace. CPU-only is an EXPLICIT
# reference/test build; a failed CUDA build or launch never selects it implicitly.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(IceHexTrace LANGUAGES CXX)
    include(CTest)
endif()
option(ICE_HEX_ENABLE_CUDA "Build the CUDA Driver-API hex tracer" ON)
option(ICE_HEX_BUILD_TESTS "Build hex tracer host and optional CUDA regression tests" ${BUILD_TESTING})
set(ICE_HEX_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(ICE_REPOSITORY_ROOT "${ICE_HEX_ROOT}/.." ABSOLUTE)

function(ice_hex_configure target)
    target_compile_features(${target} PUBLIC cxx_std_20)
    set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
    target_include_directories(${target} PUBLIC "${ICE_HEX_ROOT}/include")
    if(WIN32)
        target_compile_definitions(${target} PUBLIC NOMINMAX WIN32_LEAN_AND_MEAN)
    endif()
    if(MSVC)
        target_compile_options(${target} PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:/utf-8;/W4;/Zc:__cplusplus;/fp:precise>")
    else()
        target_compile_options(${target} PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-fno-fast-math;-ffp-contract=off>")
    endif()
endfunction()

add_library(ice_hex_trace_core STATIC src/hex_trace.cpp)
ice_hex_configure(ice_hex_trace_core)
add_executable(ice_trace_hex_cpu apps/trace_hex.cpp)
target_link_libraries(ice_trace_hex_cpu PRIVATE ice_hex_trace_core)
ice_hex_configure(ice_trace_hex_cpu)

if(ICE_HEX_ENABLE_CUDA)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 86 CACHE STRING "CUDA architectures")
    endif()
    enable_language(CUDA)
    find_package(CUDAToolkit 13.0 REQUIRED)
    # Reuse the existing Driver-API owners without changing their source.
    if(TARGET rainbow_drop_trace_host)
        set(ice_hex_driver_target rainbow_drop_trace_host)
    else()
        add_library(ice_hex_driver STATIC
            "${ICE_REPOSITORY_ROOT}/src/cuda_driver.cpp"
            "${ICE_REPOSITORY_ROOT}/src/cuda_module.cpp"
            "${ICE_REPOSITORY_ROOT}/src/cuda_error.cpp"
            "${ICE_REPOSITORY_ROOT}/src/read_binary_file.cpp")
        target_include_directories(ice_hex_driver PUBLIC "${ICE_REPOSITORY_ROOT}/include")
        target_link_libraries(ice_hex_driver PUBLIC CUDA::cuda_driver)
        ice_hex_configure(ice_hex_driver)
        set(ice_hex_driver_target ice_hex_driver)
    endif()
    add_library(ice_hex_trace_cuda STATIC src/hex_trace_cuda.cpp)
    target_link_libraries(ice_hex_trace_cuda PUBLIC ice_hex_trace_core ${ice_hex_driver_target})
    ice_hex_configure(ice_hex_trace_cuda)
    add_library(ice_hex_trace_module OBJECT device/hex_trace.cu)
    target_include_directories(ice_hex_trace_module PRIVATE "${ICE_HEX_ROOT}/include")
    set_target_properties(ice_hex_trace_module PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON
        CUDA_EXTENSIONS OFF CUDA_FATBIN_COMPILATION ON CUDA_SEPARABLE_COMPILATION OFF CUDA_RUNTIME_LIBRARY None)
    target_compile_options(ice_hex_trace_module PRIVATE
        "$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--ftz=false;--prec-div=true;--prec-sqrt=true;--fmad=false>")
    if(MSVC)
        target_compile_options(ice_hex_trace_module PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/utf-8>)
    endif()
    function(ice_hex_stage executable)
        add_custom_target(${executable}_stage
            COMMAND "${CMAKE_COMMAND}" -E make_directory "$<TARGET_FILE_DIR:${executable}>/modules"
            COMMAND "${CMAKE_COMMAND}" -E copy_if_different "$<TARGET_OBJECTS:ice_hex_trace_module>"
                "$<TARGET_FILE_DIR:${executable}>/modules/hex_trace.fatbin"
            COMMAND_EXPAND_LISTS VERBATIM)
        add_dependencies(${executable}_stage ice_hex_trace_module)
        add_dependencies(${executable} ${executable}_stage)
    endfunction()
    add_executable(ice_trace_hex apps/trace_hex.cpp)
    target_compile_definitions(ice_trace_hex PRIVATE ICE_HEX_HAS_CUDA=1)
    target_link_libraries(ice_trace_hex PRIVATE ice_hex_trace_cuda)
    ice_hex_configure(ice_trace_hex)
    ice_hex_stage(ice_trace_hex)
endif()

if(ICE_HEX_BUILD_TESTS)
    enable_testing()
    add_executable(ice_hex_host_tests tests/hex_trace_tests.cpp)
    target_link_libraries(ice_hex_host_tests PRIVATE ice_hex_trace_core)
    ice_hex_configure(ice_hex_host_tests)
    add_test(NAME ice_hex_host COMMAND ice_hex_host_tests)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(NAME ice_hex_cli COMMAND "${Python3_EXECUTABLE}"
            "${ICE_HEX_ROOT}/tests/cli_test.py" "$<TARGET_FILE:ice_trace_hex_cpu>")
    endif()
    if(ICE_HEX_ENABLE_CUDA)
        add_executable(ice_hex_cuda_tests tests/hex_trace_cuda_tests.cpp)
        target_link_libraries(ice_hex_cuda_tests PRIVATE ice_hex_trace_cuda)
        ice_hex_configure(ice_hex_cuda_tests)
        ice_hex_stage(ice_hex_cuda_tests)
        add_test(NAME ice_hex_cuda COMMAND ice_hex_cuda_tests)
        # No "success" or skip when CUDA is missing. Tests are opt-in at build time.
    endif()
endif()
