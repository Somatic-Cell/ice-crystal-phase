@echo off
setlocal
pushd "%~dp0"
if errorlevel 1 exit /b 1
if exist ".venv\Scripts\python.exe" (
    ".venv\Scripts\python.exe" "tools\generate_phase_dataset.py" %*
) else (
    python "tools\generate_phase_dataset.py" %*
)
set "RESULT=%ERRORLEVEL%"
popd
exit /b %RESULT%
