# ルートの if(BUILD_TESTING) 内の add_subdirectory(tests) から読む．
# 既存の Rainbow / rainbow_gpu_smoke はそのまま残す．

# 実装を複製せず，同じ CUDA 基盤を二つのテスト executable から参照する．
set(cuda_test_sources
    "${PROJECT_SOURCE_DIR}/src/cuda_driver.cpp"
    "${PROJECT_SOURCE_DIR}/src/cuda_module.cpp"
    "${PROJECT_SOURCE_DIR}/src/cuda_error.cpp"
    "${PROJECT_SOURCE_DIR}/src/read_binary_file.cpp"
)

# field と phase のビルド構造は同一なので，テスト名だけを差し替えて登録する．
foreach(test_suite IN ITEMS field phase)
    add_executable(rainbow_${test_suite}_cpu_tests
        test_${test_suite}_cpu.cpp
    )
    target_include_directories(rainbow_${test_suite}_cpu_tests PRIVATE
        "${PROJECT_SOURCE_DIR}/include"
    )

    add_library(rainbow_${test_suite}_test_module OBJECT
        ${test_suite}_test_kernel.cu
    )
    set_target_properties(rainbow_${test_suite}_test_module PROPERTIES
        CUDA_FATBIN_COMPILATION ON
    )
    target_include_directories(rainbow_${test_suite}_test_module PRIVATE
        "${PROJECT_SOURCE_DIR}/include"
    )
    # 非正規化数を使う phase の境界検査のため，FTZ を無効にする．
    target_compile_options(rainbow_${test_suite}_test_module PRIVATE
        $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--ftz=false>
    )

    add_executable(rainbow_${test_suite}_cuda_tests
        test_${test_suite}_cuda.cpp
        ${cuda_test_sources}
    )
    target_include_directories(rainbow_${test_suite}_cuda_tests PRIVATE
        "${PROJECT_SOURCE_DIR}/include"
    )
    target_link_libraries(rainbow_${test_suite}_cuda_tests PRIVATE CUDA::cuda_driver)
    add_dependencies(rainbow_${test_suite}_cuda_tests rainbow_${test_suite}_test_module)

    foreach(test_target IN ITEMS
            rainbow_${test_suite}_cpu_tests rainbow_${test_suite}_cuda_tests)
        if(MSVC)
            target_compile_options(${test_target} PRIVATE
                /utf-8 /Zc:__cplusplus /W4 /fp:precise)
            target_compile_definitions(${test_target} PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            # TwoSum の引き算を代数的に消去・再結合させない．
            target_compile_options(${test_target} PRIVATE -fno-fast-math -ffp-contract=off)
        endif()
    endforeach()
    if(MSVC)
        target_compile_options(rainbow_${test_suite}_test_module PRIVATE
            $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:-Xcompiler=/utf-8>
            $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:-Xcompiler=/fp:precise>
        )
    endif()

    add_test(NAME rainbow_${test_suite}_cpu
        COMMAND rainbow_${test_suite}_cpu_tests)
    add_test(NAME rainbow_${test_suite}_cuda
        COMMAND rainbow_${test_suite}_cuda_tests
            "$<TARGET_OBJECTS:rainbow_${test_suite}_test_module>"
        COMMAND_EXPAND_LISTS)
    set_tests_properties(rainbow_${test_suite}_cpu PROPERTIES LABELS "cpu;numerics")
    set_tests_properties(rainbow_${test_suite}_cuda PROPERTIES LABELS "gpu;cuda;numerics")
endforeach()

include("${CMAKE_CURRENT_LIST_DIR}/raindrop_trace_tests.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/patch_build_tests.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/patch_query_tests.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/optical_interpolation_tests.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/patch_optics_tests.cmake")