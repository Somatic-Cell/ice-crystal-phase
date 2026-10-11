@REM  .\build\hex_phase\Release\ice_generate_phase.exe ^
@REM      --out out\phase_m2_fine ^
@REM      --orientations hex_phase\configs\two_orientations.csv ^
@REM      --radius-mm 0.1 ^
@REM      --length-mm 0.2 ^
@REM      --ki 1 -0.4 0.3 ^
@REM      --ior 1.31 ^
@REM      --wavelength-nm 550 ^
@REM      --theta 72 ^
@REM      --phi 144 ^
@REM      --cdf-theta 72 ^
@REM      --cdf-phi 144

.\build\hex_atlas\Release\ice_precompute_atlas.exe ^
    --out out\atlas_550 ^
    --radius-mm 0.1 ^
    --length-mm 0.2 ^
    --ior 1.31 ^
    --wavelength-nm 550 ^
    --incident-theta-intervals 6 ^
    --incident-phi-intervals 3 ^
    --theta 900 ^
    --phi 1800 ^
    --area-subdivisions 2

python hex_atlas\python\compose_atlas.py ^
    --atlas out\atlas_550 ^
    --config hex_atlas\configs\horizontal_gaussian.json ^
    --out out\deterministic_phase