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

if not exist outputs\wave_optics mkdir outputs\wave_optics

build\Debug\rainbow_trace.exe ^
    build\Debug\modules\raindrop_trace.optixir ^
    outputs\wave_optics\drop_1mm_x.csv ^
    --sphere ^
    --radius-mm 0.4 ^
    --grid 513 ^
    --inclination-deg 20 ^
    --polarization x ^
    --patch-module build\Debug\modules\patch_build.fatbin ^
    --patch-csv outputs\wave_optics\drop_1mm_x_patches.csv ^
    --query-module build\Debug\modules\patch_query.optixir ^
    --query-theta 360 ^
    --query-phi 720 ^
    --query-csv outputs\wave_optics\drop_1mm_x_queries.csv ^
    --query-hits-csv outputs\wave_optics\drop_1mm_x_hits.csv ^
    --optics-module build\Debug\modules\patch_optics.fatbin ^
    --optics-csv outputs\wave_optics\drop_1mm_x_optics.csv ^
    --wave-csv outputs\wave_optics\drop_1mm_x_wave.csv ^
    --focal-offsets 0,0,0,0 > outputs\wave_optics\run.log 2>&1

echo Exit code: %ERRORLEVEL%