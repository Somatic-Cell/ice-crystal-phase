cmake -S . -B build -DBUILD_TESTING=ON -DICE_BUILD_HEX_TRACE=ON -DICE_HEX_ENABLE_CUDA=ON -DICE_BUILD_HEX_PHASE=ON -DICE_PHASE_ENABLE_CUDA=ON -DICE_PHASE_BUILD_TESTS=ON

cmake --build build --config Release --target ice_generate_phase ice_generate_phase_cpu ice_phase_host_tests ice_phase_cuda_tests

ctest --test-dir build -C Release -R "^ice_phase_" --output-on-failure