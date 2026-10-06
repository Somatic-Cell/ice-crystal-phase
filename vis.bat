for %S in (incoherent path focal diffraction) do (
    .venv\Scripts\python.exe tools\plot_wave_optics.py ^
        outputs\patch_audit_a1_g129\wave.csv ^
        --stage %S --component total ^
        --vmin 1e-5 --vmax 10 ^
        --out outputs\wave_optics\%S.png
)