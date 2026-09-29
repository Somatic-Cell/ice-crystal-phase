# root: cmake/patch_build.cmake の直後から読み込む．既存の host library に追加する．
target_sources(rainbow_drop_trace_host PRIVATE
    "${PROJECT_SOURCE_DIR}/src/patch_query.cpp"
)
add_dependencies(rainbow_trace rainbow_patch_query_module)

# 明示的な fmaf/fma は有効．途中の式の再結合・暗黙の FMA 化を避ける．
target_compile_options(rainbow_patch_query_module PRIVATE
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--ftz=false>
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--prec-div=true>
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--prec-sqrt=true>
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--fmad=false>
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--generate-line-info>
)
