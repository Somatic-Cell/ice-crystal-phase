rmdir /S build
mkdir build
cd build
cmake .. -DBUILD_TESTING=ON
cmake --build . --config Debug --verbose --target rainbow_trace
ctest --test-dir . -C Debug --output-on-failure --no-tests=error
cd ..