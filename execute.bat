@REM .\build\Debug\Rainbow.exe

if not exist outputs mkdir outputs

@REM build\Debug\rainbow_trace.exe ^
@REM     build\Debug\modules\raindrop_trace.optixir ^
@REM     outputs\drop_1mm_x.csv ^
@REM     --radius-mm 1.0 ^
@REM     --grid 129 ^
@REM     --inclination-deg 20 ^
@REM     --polarization x ^
@REM     --patch-module build\Debug\modules\patch_build.fatbin ^
@REM     --patch-csv outputs\drop_1mm_x_patches.csv ^
@REM     --query-module build\Debug\modules\patch_query.optixir ^
@REM     --query-theta 90 ^
@REM     --query-phi 180 ^
@REM     --query-csv outputs\drop_1mm_x_queries.csv ^
@REM     --query-hits-csv outputs\drop_1mm_x_hits.csv ^
@REM     > outputs\drop_diagnostics.log 2>&1

@REM if not exist outputs\aspherical mkdir outputs\aspherical

@REM build\Debug\rainbow_trace.exe ^
@REM     build\Debug\modules\raindrop_trace.optixir ^
@REM     outputs\aspherical\drop_1mm_y.csv ^
@REM     --radius-mm 1.0 ^
@REM     --grid 513 ^
@REM     --inclination-deg 20 ^
@REM     --polarization x ^
@REM     --patch-module build\Debug\modules\patch_build.fatbin ^
@REM     --patch-csv outputs\aspherical\drop_1mm_y_patches.csv ^
@REM     --query-module build\Debug\modules\patch_query.optixir ^
@REM     --query-theta 360 ^
@REM     --query-phi 720 ^
@REM     --query-csv outputs\aspherical\drop_1mm_y_queries.csv ^
@REM     --query-hits-csv outputs\aspherical\drop_1mm_y_hits.csv ^
@REM     --optics-module build\Debug\modules\patch_optics.fatbin ^
@REM     --optics-csv outputs\aspherical\drop_1mm_y_optics.csv ^
@REM     --wave-csv outputs\aspherical\drop_1mm_y_wave.csv ^
@REM     --focal-offsets 0,0,0,0 > outputs\aspherical\run.log 2>&1

@REM echo Exit code: %ERRORLEVEL%

@REM if not exist outputs\patch_audit_a1_g129 mkdir outputs\patch_audit_a1_g129

@REM build\Debug\rainbow_trace.exe ^
@REM     build\Debug\modules\raindrop_trace.optixir ^
@REM     outputs\patch_audit_a1_g129\vertices.csv ^
@REM     --radius-mm 1.0 ^
@REM     --wavelength-nm 700 ^
@REM     --ior 1.3314 ^
@REM     --grid 3001 ^
@REM     --inclination-deg 20 ^
@REM     --polarization x ^
@REM     --patch-module build\Debug\modules\patch_build.fatbin ^
@REM     --query-module build\Debug\modules\patch_query.optixir ^
@REM     --query-theta 360 ^
@REM     --query-phi 720 ^
@REM     --query-csv outputs\patch_audit_a1_g129\queries.csv ^
@REM     --query-hits-csv outputs\patch_audit_a1_g129\hits.csv ^
@REM     --optics-module build\Debug\modules\patch_optics.fatbin ^
@REM     --optics-csv outputs\patch_audit_a1_g129\optics.csv ^
@REM     --wave-csv outputs\patch_audit_a1_g129\wave.csv ^
@REM     --focal-offsets 0,0,0,0 ^
@REM     > outputs\patch_audit_a1_g129\run.log 2>&1

@REM echo Exit code: %ERRORLEVEL%

@REM .venv\Scripts\python.exe tools\analyze_patch_failures.py ^
@REM     outputs\patch_audit_a1_g129\patch_failure_report.json ^
@REM     --out outputs\patch_audit_a1_g129\exact_audit

@REM echo Audit exit: %ERRORLEVEL%

if not exist outputs\unpolarized_a1_g129 mkdir outputs\unpolarized_a1_g129

build\Release\rainbow_generate.exe ^
    --out datasets\drop_a1_i20_700nm_q1800x3600_c90x1800 ^
    --radius-mm 1.0 ^
    --wavelength-nm 700 ^
    --temperature-c 20 ^
    --pressure-pa 101325 ^
    --grid 3001 ^
    --inclination-deg 20 ^
    --query-theta 1800 ^
    --query-phi 3600 ^
    --cdf-theta 900 ^
    --cdf-phi 1800 ^
    --stage diffraction ^
    --focal-offsets 0,0,0,0 ^
    --allow-underresolved

echo Exit code: %ERRORLEVEL%
