#pragma once

namespace rainbow
{
// Host-only material evaluation. All wavelengths are vacuum wavelengths in nm.
// This material service is independent of CUDA, OptiX, ray tracing and CDF size.
struct WaterOpticalState
{
    double temperature_kelvin = 293.15;
    double pressure_pascal = 101325.0;
    double density_kg_m3 = 0.0;
    double density_pressure_residual_pascal = 0.0;
    double density_bracket_width_kg_m3 = 0.0;
    double refractive_index = 0.0; // real absolute phase index; not group index
};

// R9-97, Eq. (1), all coefficients of Table 1. Does not model absorption.
// Endorsed input range: 200..1100 nm, 261.15..773.15 K, 0..1060 kg/m^3.
[[nodiscard]] double water_refractive_index_r9_97(
    double vacuum_wavelength_nm, double temperature_kelvin, double density_kg_m3);

// IAPWS-95 pressure EOS -> liquid density -> R9-97 real index.
// Intentionally bounded application domain: 0..60 C and 50 kPa..100 MPa.
// Not a general water/steam phase-selection or saturation solver. Rejects inputs
// outside this contract instead of extrapolating or substituting a density fit.
[[nodiscard]] WaterOpticalState evaluate_water_optics(
    double vacuum_wavelength_nm,
    double temperature_celsius = 20.0,
    double pressure_pascal = 101325.0);
} // namespace rainbow
