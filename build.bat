rmdir /S build
mkdir build
cd build
cmake .. -DBUILD_TESTING=ON
cmake --build . --config Debug --verbose
ctest --test-dir build -C Debug --output-on-failure --no-tests=error
cd ..