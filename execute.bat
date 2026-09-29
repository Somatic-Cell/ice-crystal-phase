@REM .\build\Debug\Rainbow.exe

if not exist outputs mkdir outputs

build\Debug\rainbow_trace.exe ^
    build\Debug\modules\raindrop_trace.optixir ^
    outputs\drop_1mm_x.csv ^
    --radius-mm 1.0 ^
    --grid 129 ^
    --inclination-deg 20 ^
    --polarization x ^
    --patch-module build\Debug\modules\patch_build.fatbin ^
    --patch-csv outputs\drop_1mm_x_patches.csv ^
    --query-module build\Debug\modules\patch_query.optixir ^
    --query-theta 90 ^
    --query-phi 180 ^
    --query-csv outputs\drop_1mm_x_queries.csv ^
    --query-hits-csv outputs\drop_1mm_x_hits.csv
    > outputs\drop_diagnostics.log 2>&1