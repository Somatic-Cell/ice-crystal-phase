@REM  .venv\Scripts\python.exe tools\plot_phase_cdf.py ^
@REM      out\deterministic_phase ^
@REM      --device auto ^
@REM      --out out\cdf_pdf_log.png ^
@REM      --report out\cdf_pdf_log.json

for /L %%i in (0,1,6)  do python tools\phase_analysis.py inspect "out\deterministic_phase\records\i000%%i\w0000" ^
        --out "analysis\m3_all\i000%%i\w0000" ^
        --label "M3: incident index %%i, wavelength 550 nm"