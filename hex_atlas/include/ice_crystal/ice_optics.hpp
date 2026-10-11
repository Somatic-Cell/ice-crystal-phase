#pragma once
#include <ostream>

namespace iceCrystal::optics
{
// Host-side monochromatic state; no wavelength lookup in the CUDA ray loop.
// n is used by the existing isotropic, nonabsorbing tracer. kappa is evaluated
// and reported but is NOT applied to that tracer. No temperature extrapolation.
struct IceState
{
    double wavelength_nm = 0;
    double n = 0;
    double kappa = 0;
    double lower_nm = 0, upper_nm = 0, log_wavelength_fraction = 0;
};
[[nodiscard]] IceState evaluate_ice_wb2008(double vacuum_wavelength_nm);
// Writes one JSON object. The caller owns the enclosing document and locale.
void write_material_json(std::ostream& out, const IceState& state,
                         double exterior_index, bool ice_material);
} // namespace iceCrystal::optics
