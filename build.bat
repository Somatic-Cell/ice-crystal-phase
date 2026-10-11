cmake -S . -B build -DICE_BUILD_HEX_ATLAS=ON -DICE_ATLAS_ENABLE_CUDA=ON -DICE_ATLAS_BUILD_TESTS=ON
cmake --build build --config Release --target ice_precompute_atlas ice_precompute_atlas_cpu ice_atlas_host_tests ice_atlas_cuda_tests
ctest --test-dir build -C Release -R "^ice_atlas_" --output-on-failure
