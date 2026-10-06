@REM rmdir /S build
@REM mkdir build
@REM cd build
@REM cmake -S .. -B . -DBUILD_TESTING=ON && ^
@REM cmake --build . --config Debug && ^
@REM ctest --test-dir . -C Debug --output-on-failure --no-tests=error
@REM cd ..

cmake -S . -B build -DBUILD_TESTING=ON && ^
cmake --build build --config Debug && ^
ctest --test-dir build -C Debug --output-on-failure --no-tests=error

.venv\Scripts\python.exe -m unittest discover ^
    -s tests -p test_patch_failure_analysis.py -v