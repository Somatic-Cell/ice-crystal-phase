# root CMakeLists.txt の 3 個の CUDA/OptiX module target 定義の後から include する．
# 旧 raygen-only smoke を本体から分離．既存 M2 はこの module に切り替える．
add_library(rainbow_optix_smoke_module OBJECT
    "${PROJECT_SOURCE_DIR}/tests/optix_smoke.cu"
)
set_target_properties(rainbow_optix_smoke_module PROPERTIES CUDA_OPTIX_COMPILATION ON)
target_include_directories(rainbow_optix_smoke_module PRIVATE
    "${PROJECT_SOURCE_DIR}/include" "${OPTIX91_INCLUDE_DIR}"
)
if(MSVC)
    target_compile_options(rainbow_optix_smoke_module PRIVATE
        $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:-Xcompiler=/utf-8>
    )
endif()

# Interval bounds / FloatPair は各演算の丸めに依存する．global fast-math を使わない．
# fmaf() で明示した FMA は有効なまま，暗黙の式の融合だけを禁止する．
target_compile_options(rainbow_raindrop_trace_module PRIVATE
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--ftz=false>
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--prec-div=true>
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--prec-sqrt=true>
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--fmad=false>
)

# function table の定義と初期化は共通 TU へ移す．旧 src/optix.cpp の定義 include は削除．
target_sources(Rainbow PRIVATE "${PROJECT_SOURCE_DIR}/src/optix_context.cpp")

# 実際の tracing executable と結合検査が同じホスト実装をリンクする．
# CUDA 初期化/エラー処理/file reader を複製しない．
add_library(rainbow_drop_trace_host STATIC
    "${PROJECT_SOURCE_DIR}/src/raindrop_tracer.cpp"
    "${PROJECT_SOURCE_DIR}/src/optix_context.cpp"
    "${PROJECT_SOURCE_DIR}/src/cuda_driver.cpp"
    "${PROJECT_SOURCE_DIR}/src/cuda_module.cpp"
    "${PROJECT_SOURCE_DIR}/src/cuda_error.cpp"
    "${PROJECT_SOURCE_DIR}/src/optix_error.cpp"
    "${PROJECT_SOURCE_DIR}/src/read_binary_file.cpp"
)

if(MSVC)
    # デバッグ情報のコピー元がディレクトリだけにならないよう，
    # コンパイラが生成する PDB の名前と配置先を明示する．
    set_target_properties(
        rainbow_drop_trace_host
        PROPERTIES
            COMPILE_PDB_NAME "rainbow_drop_trace_host"
            COMPILE_PDB_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/pdb"
    )
endif()


target_include_directories(rainbow_drop_trace_host PUBLIC
    "${PROJECT_SOURCE_DIR}/include" "${OPTIX91_INCLUDE_DIR}"
)
target_link_libraries(rainbow_drop_trace_host PUBLIC CUDA::cuda_driver)
if(UNIX)
    target_link_libraries(rainbow_drop_trace_host PRIVATE ${CMAKE_DL_LIBS})
endif()
if(MSVC)
    target_compile_options(rainbow_drop_trace_host PRIVATE /utf-8 /W4 /Zc:__cplusplus /fp:precise)
    target_compile_definitions(rainbow_drop_trace_host PUBLIC NOMINMAX WIN32_LEAN_AND_MEAN)
else()
    target_compile_options(rainbow_drop_trace_host PRIVATE -ffp-contract=off -fno-fast-math)
endif()

add_executable(rainbow_trace "${PROJECT_SOURCE_DIR}/apps/trace_raindrop.cpp")
target_link_libraries(rainbow_trace PRIVATE rainbow_drop_trace_host)
if(MSVC)
    target_compile_options(rainbow_trace PRIVATE /utf-8 /W4 /Zc:__cplusplus /fp:precise)
endif()
add_dependencies(rainbow_trace rainbow_raindrop_trace_module)
# CLI は IR のパスを明示的に受け取る．root の stage_modules が .exe 隣へ配置する．
