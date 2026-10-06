@echo off
setlocal
cd /d "%~dp0"
rem Explicit regression run in a separate build tree; never runs from build.bat.
cmake -S . -B build-verify -DBUILD_TESTING=ON -DRAINBOW_ENABLE_DIAGNOSTICS=ON %*
if errorlevel 1 exit /b %ERRORLEVEL%
cmake --build build-verify --config Release
if errorlevel 1 exit /b %ERRORLEVEL%
ctest --test-dir build-verify -C Release --output-on-failure --no-tests=error
exit /b %ERRORLEVEL%
