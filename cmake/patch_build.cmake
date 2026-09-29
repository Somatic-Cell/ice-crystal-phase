# Include from root, immediately AFTER cmake/raindrop_trace.cmake.
# All referenced targets have already been defined at that point.
target_sources(rainbow_drop_trace_host PRIVATE
    "${PROJECT_SOURCE_DIR}/src/patch_accel.cpp"
)
add_dependencies(rainbow_trace rainbow_patch_build_module)

# FP32 interval sign certificates must not be invalidated by reassociation/FTZ.
# Explicit fmaf calls remain fused; implicit multiplication/addition does not.
target_compile_options(rainbow_patch_build_module PRIVATE
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--ftz=false>
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--prec-div=true>
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--prec-sqrt=true>
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--fmad=false>
)
# No OptiX include directory is added to the ordinary CUDA kernel target.
# The existing stage_modules target already copies patch_build.fatbin.
