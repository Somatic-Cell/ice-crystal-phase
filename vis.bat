@REM  for %S in (incoherent path focal diffraction) do (
@REM      .venv\Scripts\python.exe tools\plot_wave_optics.py ^
@REM          outputs\unpolarized_a1_g129\wave.csv ^
@REM          --stage %S --component total ^
@REM          --vmin 1e-5 --vmax 10 ^
@REM          --out outputs\wave_optics\%S.png
@REM  )

.venv\Scripts\python.exe tools\plot_phase_cdf.py ^
    datasets\drop_a1_i20_700nm_q1800x3600_c90x1800 ^
    --device auto ^
    --out outputs\aspherical\cdf_pdf_log.png ^
    --report outputs\aspherical\cdf_pdf_log.json

@REM  .venv\Scripts\python.exe tools\plot_phase_cdf.py ^
@REM      datasets\drop_a1_i20_700nm_q900x1800_c450x900 ^
@REM      --device cuda ^
@REM      --theta-range 70 130 ^
@REM      --phi-range -25 25 ^
@REM      --log-limits -6 2 ^
@REM      --out outputs\cdf_pdf_center.png ^
@REM      --cell-image outputs\cdf_pdf_center_cells.png