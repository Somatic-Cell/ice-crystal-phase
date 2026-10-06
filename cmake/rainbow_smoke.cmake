# Included only with BUILD_TESTING=ON. The legacy app is a smoke test, not a solver.
# The unchanged main uses GetModuleFileNameW; do not claim cross-platform support.
if(NOT WIN32)
    message(STATUS "Windows-only legacy rainbow_gpu_smoke is not registered on this platform.")
    return()
endif()
add_library(rainbow_optix_smoke_module OBJECT
    "${PROJECT_SOURCE_DIR}/tests/optix_smoke.cu")
set_target_properties(rainbow_optix_smoke_module PROPERTIES
    CUDA_OPTIX_COMPILATION ON CUDA_SEPARABLE_COMPILATION OFF CUDA_RUNTIME_LIBRARY None)
target_include_directories(rainbow_optix_smoke_module PRIVATE
    "${PROJECT_SOURCE_DIR}/include" "${OPTIX91_INCLUDE_DIR}")
if(MSVC)
    target_compile_options(rainbow_optix_smoke_module PRIVATE
        $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:-Xcompiler=/utf-8>)
endif()

add_executable(rainbow_gpu_smoke_tests
    "${PROJECT_SOURCE_DIR}/tests/smoke/main.cpp"
    "${PROJECT_SOURCE_DIR}/tests/smoke/optix.cpp")
target_link_libraries(rainbow_gpu_smoke_tests PRIVATE rainbow::solver)
rainbow_configure_host(rainbow_gpu_smoke_tests)
# Separate directory: do not race the solver's stage target on the same files.
set_target_properties(rainbow_gpu_smoke_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/tests/smoke")
add_custom_target(rainbow_stage_smoke_modules
    COMMAND "${CMAKE_COMMAND}" -E make_directory
        "$<TARGET_FILE_DIR:rainbow_gpu_smoke_tests>/modules"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different
        "$<TARGET_OBJECTS:rainbow_patch_build_module>"
        "$<TARGET_FILE_DIR:rainbow_gpu_smoke_tests>/modules/patch_build.fatbin"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different
        "$<TARGET_OBJECTS:rainbow_optix_smoke_module>"
        "$<TARGET_FILE_DIR:rainbow_gpu_smoke_tests>/modules/optix_smoke.optixir"
    COMMAND_EXPAND_LISTS VERBATIM)
add_dependencies(rainbow_stage_smoke_modules rainbow_patch_build_module rainbow_optix_smoke_module)
add_dependencies(rainbow_gpu_smoke_tests rainbow_stage_smoke_modules)
add_test(NAME rainbow_gpu_smoke COMMAND rainbow_gpu_smoke_tests)
set_tests_properties(rainbow_gpu_smoke PROPERTIES LABELS "gpu;cuda;optix")
