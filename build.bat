@echo off
setlocal
cd /d "%~dp0"
rem Build only. No clean build, no tests, no simulation, no visualization.
cmake -S . -B build -DBUILD_TESTING=OFF -DRAINBOW_ENABLE_DIAGNOSTICS=OFF %*
if errorlevel 1 exit /b %ERRORLEVEL%
cmake --build build --config Release --target rainbow_generate
exit /b %ERRORLEVEL%