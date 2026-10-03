rmdir /S build
mkdir build
cd build
cmake -S .. -B . -DBUILD_TESTING=ON && ^
cmake --build . --config Debug && ^
ctest --test-dir . -C Debug --output-on-failure --no-tests=error
cd ..