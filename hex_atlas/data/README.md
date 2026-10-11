# Frozen ice optical constants

Warren, S. G., and R. E. Brandt (2008), *Optical constants of ice from the
ultraviolet to the microwave: A revised compilation*, JGR 113, D14220.
DOI: 10.1029/2007JD009744.

Author source: https://www.atmos.washington.edu/ice_optical_constants/
Table: https://www.atmos.washington.edu/ice_optical_constants/IOP_2008_ASCIItable.dat

`wb2008_199_1100nm.dat` preserves the numeric entries from 0.199 through 1.100 um
(79 rows). The first row brackets the explicitly supported lower bound, 200 nm.
The data are compiled into `src/ice_optics_table.inc`; no file/network lookup is
performed in the CUDA ray loop or at run time. Unit conversion: um * 1000 = nm.
The raw subset SHA256 is
`71235dd08b3f65fc5fdf638101577d4da12af9cd725da21fa1677c0f91bfddcb`.

Following the author's interpolation instructions: n is linear in log(wavelength),
log(kappa) is linear in log(wavelength). Knots are returned without interpolation.
No extrapolation beyond 200..1100 nm. This is the ice-Ih near-melting compilation,
not an arbitrary temperature/pressure model. The real n is used in the existing
isotropic nonabsorbing geometric-optics solver. kappa is evaluated and reported
but **not used**; the solver's kappa remains zero. Birefringence, diffraction,
absorption and inter-path interference have not been added by this patch.

This is numerical data, not copied paper text or an imported material library.
Please retain this attribution when reusing the numeric subset.
