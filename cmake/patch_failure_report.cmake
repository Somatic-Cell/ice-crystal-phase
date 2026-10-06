# Include from ROOT CMakeLists.txt, after cmake/patch_optics.cmake.
# Diagnostics are host-only. No CUDA/OptiX module or numerical option changes.
target_sources(rainbow_drop_trace_host PRIVATE
    "${PROJECT_SOURCE_DIR}/src/patch_failure_report.cpp"
    "${PROJECT_SOURCE_DIR}/src/patch_failure_capture.cpp"
)
# FP32/FP64 enclosure audit must not contract/reassociate interval expressions.
if(MSVC)
    set_property(SOURCE "${PROJECT_SOURCE_DIR}/src/patch_failure_report.cpp"
        APPEND PROPERTY COMPILE_OPTIONS /fp:strict)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set_property(SOURCE "${PROJECT_SOURCE_DIR}/src/patch_failure_report.cpp"
        APPEND PROPERTY COMPILE_OPTIONS -fno-fast-math -ffp-contract=off)
endif()
