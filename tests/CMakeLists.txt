# 親の CMakeLists.txt の if(BUILD_TESTING) 内から add_subdirectory(tests) で読む．
# 既存の Rainbow / rainbow_gpu_smoke は変更しない．

# CPU 側は純粋な C++ executable．CUDA / OptiX のライブラリをリンクしない．
add_executable(rainbow_field_cpu_tests
    test_field_cpu.cpp
)

target_include_directories(rainbow_field_cpu_tests PRIVATE
    "${PROJECT_SOURCE_DIR}/include"
)

# 通常 CUDA の検証用モジュール．本番用 patch_build.cu に検査コードを追加しない．
add_library(rainbow_field_test_module OBJECT
    field_test_kernel.cu
)

set_target_properties(rainbow_field_test_module PROPERTIES
    CUDA_FATBIN_COMPILATION ON
)

target_include_directories(rainbow_field_test_module PRIVATE
    "${PROJECT_SOURCE_DIR}/include"
)

# 既存実装のソースを参照し，独自の CUDA 初期化・module loader を複製しない．
# CudaContext がまだ cuda_driver.cpp にあるため，同ファイルをビルド対象にする．
# 同じソースが executable ごとにコンパイルされるが，実装のコピーは作らない．
add_executable(rainbow_field_cuda_tests
    test_field_cuda.cpp
    "${PROJECT_SOURCE_DIR}/src/cuda_driver.cpp"
    "${PROJECT_SOURCE_DIR}/src/cuda_module.cpp"
    "${PROJECT_SOURCE_DIR}/src/cuda_error.cpp"
    "${PROJECT_SOURCE_DIR}/src/read_binary_file.cpp"
)

target_include_directories(rainbow_field_cuda_tests PRIVATE
    "${PROJECT_SOURCE_DIR}/include"
)

target_link_libraries(rainbow_field_cuda_tests PRIVATE
    CUDA::cuda_driver
)

# 警告・文字コード設定は新しい target にも必要．root の PRIVATE 設定は継承されない．
if(MSVC)
    foreach(test_target IN ITEMS rainbow_field_cpu_tests rainbow_field_cuda_tests)
        target_compile_options(${test_target} PRIVATE
            /utf-8
            /Zc:__cplusplus
            /W4
        )
        target_compile_definitions(${test_target} PRIVATE
            NOMINMAX
            WIN32_LEAN_AND_MEAN
        )
    endforeach()
    target_compile_options(rainbow_field_test_module PRIVATE
        $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:-Xcompiler=/utf-8>
    )
endif()

add_dependencies(rainbow_field_cuda_tests rainbow_field_test_module)

add_test(NAME rainbow_field_cpu
    COMMAND rainbow_field_cpu_tests
)

# TARGET_OBJECTS は1個の .cu が生成する1個の .fatbin のパスへ展開する．
# CTest から生成元を直接指定するので，テスト用モジュールのコピー処理は不要．
add_test(NAME rainbow_field_cuda
    COMMAND rainbow_field_cuda_tests "$<TARGET_OBJECTS:rainbow_field_test_module>"
    COMMAND_EXPAND_LISTS
)

set_tests_properties(rainbow_field_cpu PROPERTIES LABELS "cpu;numerics")
set_tests_properties(rainbow_field_cuda PROPERTIES LABELS "gpu;cuda;numerics")
