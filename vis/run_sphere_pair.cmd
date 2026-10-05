@echo off
setlocal
rem Run from any current directory. This script lives in <repository>\vis\.
pushd "%~dp0.." || exit /b 1
set "EXE=build\Debug\rainbow_trace.exe"
set "MODULES=build\Debug\modules"
set "OUT=outputs\sphere_validation"
rem Initial comparison resolution, NOT a convergence guarantee.
set "GRID=513"
set "THETA=1800"
set "PHI=64"
if not exist "%EXE%" (
    echo Missing %EXE%. Build rainbow_trace first.
    popd
    exit /b 1
)
if not exist "%OUT%" mkdir "%OUT%"
if not exist "%OUT%" (
    echo Cannot create %OUT%.
    popd
    exit /b 1
)
for %%P in (x y) do (
    echo Running sphere with incident polarization %%P ...
    "%EXE%" ^
        "%MODULES%\raindrop_trace.optixir" ^
        "%OUT%\sphere_%%P_vertices.csv" ^
        --sphere --radius-mm 0.4 --wavelength-nm 700 --ior 1.3314 ^
        --grid %GRID% --inclination-deg 0 --polarization %%P ^
        --patch-module "%MODULES%\patch_build.fatbin" ^
        --query-module "%MODULES%\patch_query.optixir" ^
        --query-theta %THETA% --query-phi %PHI% ^
        --query-csv "%OUT%\sphere_%%P_queries.csv" ^
        --optics-module "%MODULES%\patch_optics.fatbin" ^
        --optics-csv "%OUT%\sphere_%%P_optics.csv" > "%OUT%\sphere_%%P.log" 2>&1
    if errorlevel 1 (
        echo Simulation failed. Inspect %OUT%\sphere_%%P.log
        popd
        exit /b 1
    )
)
echo Both sphere runs finished. This is not a physical-validation certificate.
popd
exit /b 0
