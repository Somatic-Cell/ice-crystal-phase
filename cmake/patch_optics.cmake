# Root: include after cmake/patch_query.cmake. Numerical device code has no
# OptiX dependency; it consumes completed query buffers using CUDA Driver API.
add_library(rainbow_patch_optics_module OBJECT
    "${PROJECT_SOURCE_DIR}/device/patch_optics.cu"
)
set_target_properties(rainbow_patch_optics_module PROPERTIES
    CUDA_FATBIN_COMPILATION ON
    CUDA_SEPARABLE_COMPILATION OFF
)
target_include_directories(rainbow_patch_optics_module PRIVATE
    "${PROJECT_SOURCE_DIR}/include"
)
target_compile_options(rainbow_patch_optics_module PRIVATE
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--ftz=false>
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--prec-div=true>
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--prec-sqrt=true>
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--fmad=false>
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--generate-line-info>
)
if(MSVC)
    target_compile_options(rainbow_patch_optics_module PRIVATE
        $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:-Xcompiler=/utf-8>
        $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:-Xcompiler=/fp:precise>
    )
endif()
target_sources(rainbow_drop_trace_host PRIVATE
    "${PROJECT_SOURCE_DIR}/src/patch_optics.cpp"
)
add_dependencies(rainbow_trace rainbow_patch_optics_module)
