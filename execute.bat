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

cmake -E env OPTIX_CACHE_MAXSIZE=0 ^
    build\Debug\rainbow_trace.exe ^
    build\Debug\modules\raindrop_trace.optixir ^
    outputs\validation_off_probe\drop_1mm_x.csv ^
    --radius-mm 1.0 ^
    --grid 129 ^
    --inclination-deg 20 ^
    --polarization x ^
    --patch-module build\Debug\modules\patch_build.fatbin ^
    --patch-csv outputs\validation_off_probe\drop_1mm_x_patches.csv ^
    --query-module build\Debug\modules\patch_query.optixir ^
    --query-theta 90 ^
    --query-phi 180 ^
    --query-csv outputs\validation_off_probe\drop_1mm_x_queries.csv ^
    --query-hits-csv outputs\validation_off_probe\drop_1mm_x_hits.csv > outputs\validation_off_probe\run.log 2>&1

echo Exit code: %ERRORLEVEL%