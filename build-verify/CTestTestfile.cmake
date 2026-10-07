# CMake generated Testfile for 
# Source directory: C:/Users/sy415/workspace/rainbow
# Build directory: C:/Users/sy415/workspace/rainbow/build-verify
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
  add_test([=[rainbow_gpu_smoke]=] "C:/Users/sy415/workspace/rainbow/build-verify/tests/smoke/Debug/rainbow_gpu_smoke_tests.exe")
  set_tests_properties([=[rainbow_gpu_smoke]=] PROPERTIES  LABELS "gpu;cuda;optix" _BACKTRACE_TRIPLES "C:/Users/sy415/workspace/rainbow/cmake/rainbow_smoke.cmake;38;add_test;C:/Users/sy415/workspace/rainbow/cmake/rainbow_smoke.cmake;0;;C:/Users/sy415/workspace/rainbow/CMakeLists.txt;148;include;C:/Users/sy415/workspace/rainbow/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
  add_test([=[rainbow_gpu_smoke]=] "C:/Users/sy415/workspace/rainbow/build-verify/tests/smoke/Release/rainbow_gpu_smoke_tests.exe")
  set_tests_properties([=[rainbow_gpu_smoke]=] PROPERTIES  LABELS "gpu;cuda;optix" _BACKTRACE_TRIPLES "C:/Users/sy415/workspace/rainbow/cmake/rainbow_smoke.cmake;38;add_test;C:/Users/sy415/workspace/rainbow/cmake/rainbow_smoke.cmake;0;;C:/Users/sy415/workspace/rainbow/CMakeLists.txt;148;include;C:/Users/sy415/workspace/rainbow/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
  add_test([=[rainbow_gpu_smoke]=] "C:/Users/sy415/workspace/rainbow/build-verify/tests/smoke/MinSizeRel/rainbow_gpu_smoke_tests.exe")
  set_tests_properties([=[rainbow_gpu_smoke]=] PROPERTIES  LABELS "gpu;cuda;optix" _BACKTRACE_TRIPLES "C:/Users/sy415/workspace/rainbow/cmake/rainbow_smoke.cmake;38;add_test;C:/Users/sy415/workspace/rainbow/cmake/rainbow_smoke.cmake;0;;C:/Users/sy415/workspace/rainbow/CMakeLists.txt;148;include;C:/Users/sy415/workspace/rainbow/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
  add_test([=[rainbow_gpu_smoke]=] "C:/Users/sy415/workspace/rainbow/build-verify/tests/smoke/RelWithDebInfo/rainbow_gpu_smoke_tests.exe")
  set_tests_properties([=[rainbow_gpu_smoke]=] PROPERTIES  LABELS "gpu;cuda;optix" _BACKTRACE_TRIPLES "C:/Users/sy415/workspace/rainbow/cmake/rainbow_smoke.cmake;38;add_test;C:/Users/sy415/workspace/rainbow/cmake/rainbow_smoke.cmake;0;;C:/Users/sy415/workspace/rainbow/CMakeLists.txt;148;include;C:/Users/sy415/workspace/rainbow/CMakeLists.txt;0;")
else()
  add_test([=[rainbow_gpu_smoke]=] NOT_AVAILABLE)
endif()
subdirs("tests")
