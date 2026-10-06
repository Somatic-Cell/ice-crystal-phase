# Only the four production device modules. Do not link their OBJECT outputs.
include_guard(GLOBAL)

function(rainbow_add_module target source kind)
    add_library(${target} OBJECT "${PROJECT_SOURCE_DIR}/${source}")
    set_target_properties(${target} PROPERTIES
        CUDA_SEPARABLE_COMPILATION OFF
        CUDA_RUNTIME_LIBRARY None)
    if(kind STREQUAL "OPTIX")
        set_target_properties(${target} PROPERTIES CUDA_OPTIX_COMPILATION ON)
        target_include_directories(${target} PRIVATE
            "${PROJECT_SOURCE_DIR}/include"
            "${PROJECT_SOURCE_DIR}/shaders"
            "${OPTIX91_INCLUDE_DIR}")
    elseif(kind STREQUAL "FATBIN")
        set_target_properties(${target} PROPERTIES CUDA_FATBIN_COMPILATION ON)
        target_include_directories(${target} PRIVATE
            "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/device")
    else()
        message(FATAL_ERROR "Unknown device-module kind: ${kind}")
    endif()
    # Preserve the baseline's interval/FloatPair arithmetic contract.
    # Explicit fma/fmaf calls are not disabled by --fmad=false.
    target_compile_options(${target} PRIVATE
        $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--ftz=false>
        $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--prec-div=true>
        $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--prec-sqrt=true>
        $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--fmad=false>)
    if(MSVC)
        target_compile_options(${target} PRIVATE
            $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:-Xcompiler=/utf-8>)
    endif()
endfunction()

rainbow_add_module(rainbow_raindrop_trace_module shaders/entry/raindrop_trace.cu OPTIX)
rainbow_add_module(rainbow_patch_build_module device/patch_build.cu FATBIN)
rainbow_add_module(rainbow_patch_query_module shaders/entry/patch_query.cu OPTIX)
rainbow_add_module(rainbow_patch_optics_module device/patch_optics.cu FATBIN)

# These two modules had line information in the baseline; retain it.
foreach(target IN ITEMS rainbow_patch_query_module rainbow_patch_optics_module)
    target_compile_options(${target} PRIVATE
        $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--generate-line-info>)
endforeach()
if(MSVC)
    target_compile_options(rainbow_patch_optics_module PRIVATE
        $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:-Xcompiler=/fp:precise>)
endif()

function(rainbow_stage_solver_modules executable stage_target)
    # A separate target runs even when a module, but not the .exe, was rebuilt.
    # TARGET_FILE_DIR does not create a reverse dependency (CMP0112 NEW).
    # No TARGET_FILE generator expression is used here: avoid a dependency cycle.
    add_custom_target(${stage_target}
        COMMAND "${CMAKE_COMMAND}" -E make_directory
            "$<TARGET_FILE_DIR:${executable}>/modules"
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different
            "$<TARGET_OBJECTS:rainbow_raindrop_trace_module>"
            "$<TARGET_FILE_DIR:${executable}>/modules/raindrop_trace.optixir"
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different
            "$<TARGET_OBJECTS:rainbow_patch_build_module>"
            "$<TARGET_FILE_DIR:${executable}>/modules/patch_build.fatbin"
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different
            "$<TARGET_OBJECTS:rainbow_patch_query_module>"
            "$<TARGET_FILE_DIR:${executable}>/modules/patch_query.optixir"
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different
            "$<TARGET_OBJECTS:rainbow_patch_optics_module>"
            "$<TARGET_FILE_DIR:${executable}>/modules/patch_optics.fatbin"
        COMMAND_EXPAND_LISTS VERBATIM)
    add_dependencies(${stage_target}
        rainbow_raindrop_trace_module rainbow_patch_build_module
        rainbow_patch_query_module rainbow_patch_optics_module)
    add_dependencies(${executable} ${stage_target})
endfunction()
