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

# 診断用：同じ patch_query.cu を PTX にもコンパイルする．
# 既存の OptiX-IR target と，本体の module 配置は変更しない．
add_library(
    rainbow_patch_query_ptx_probe
    OBJECT
        "${PROJECT_SOURCE_DIR}/shaders/entry/patch_query.cu"
)

set_target_properties(
    rainbow_patch_query_ptx_probe
    PROPERTIES
        CUDA_PTX_COMPILATION ON
)

# include・精度関連オプション・マクロ定義は既存 target と揃える．
# --use_fast_math などの別の数値条件を追加しない．
target_include_directories(
    rainbow_patch_query_ptx_probe
    PRIVATE
        "$<TARGET_PROPERTY:rainbow_patch_query_module,INCLUDE_DIRECTORIES>"
)

target_compile_options(
    rainbow_patch_query_ptx_probe
    PRIVATE
        "$<TARGET_PROPERTY:rainbow_patch_query_module,COMPILE_OPTIONS>"
)

target_compile_definitions(
    rainbow_patch_query_ptx_probe
    PRIVATE
        "$<TARGET_PROPERTY:rainbow_patch_query_module,COMPILE_DEFINITIONS>"
)

add_dependencies(
    rainbow_patch_query_optix_tests
    rainbow_patch_query_ptx_probe
)

# 新しい C++ テストは書かず，同じ executable に PTX を渡す．
add_test(
    NAME rainbow_patch_query_optix_ptx_probe
    COMMAND
        rainbow_patch_query_optix_tests
        "$<TARGET_OBJECTS:rainbow_patch_build_module>"
        "$<TARGET_OBJECTS:rainbow_patch_query_ptx_probe>"
    COMMAND_EXPAND_LISTS
)

# 今回の比較では，両方ともコンパイルキャッシュを使わない．
# キャッシュの内容を削除する設定ではない．
set_tests_properties(
    rainbow_patch_query_optix
    rainbow_patch_query_optix_ptx_probe
    PROPERTIES
        ENVIRONMENT "OPTIX_CACHE_MAXSIZE=0"
)
