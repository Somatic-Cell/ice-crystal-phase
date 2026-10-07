cmake -S . -B build -DBUILD_TESTING=ON -DICE_BUILD_HEX_TRACE=ON -DICE_HEX_ENABLE_CUDA=ON -DICE_HEX_BUILD_TESTS=ON

cmake --build build --config Release --target ice_trace_hex ice_trace_hex_cpu ice_hex_host_tests ice_hex_cuda_tests

ctest --test-dir build -C Release -R "^ice_hex_" --output-on-failure