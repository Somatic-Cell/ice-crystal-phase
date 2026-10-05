"""Full Mie *intensity* reference for precisely matching incident polarization.

miepython intensities(..., norm='qsca') returns (parallel, perpendicular).
Multiplication by pi*radius**2 converts this to differential cross section.
The four ray families in the C++ model are NOT four multipole terms: n_pole=0
is used intentionally. We do not compare complex phases across conventions.
"""
from __future__ import annotations

import math
from typing import Any
import numpy as np

from optics_io import IntensityView, metadata_float, metadata_vector


def incident_projection_weights(view: IntensityView) -> tuple[np.ndarray, np.ndarray]:
    """Fractions in (s,p_in), including arbitrary complex incident Jones states."""
    g = view.grid
    if view.is_unpolarized:
        return np.full(g.shape, 0.5), np.full(g.shape, 0.5)
    # e1 = wi cross e0. At azimuth phi:
    # s = -sin(phi)*e0 + cos(phi)*e1; p_in = wi cross s.
    j = g.jones / math.sqrt(g.incident_norm2)
    s = -np.sin(g.phi)*j[0] + np.cos(g.phi)*j[1]
    p = -np.cos(g.phi)*j[0] - np.sin(g.phi)*j[1]
    fs = np.broadcast_to(np.abs(s)**2, g.shape).copy()
    fp = np.broadcast_to(np.abs(p)**2, g.shape).copy()
    if not np.allclose(fs+fp, 1, rtol=0, atol=1e-12):
        raise ValueError("Incident polarization decomposition is not normalized")
    return fs, fp


def calculate_mie(view: IntensityView, *, backend=None) -> tuple[np.ndarray, dict[str, Any]]:
    if not all(g.is_sphere for g in view.grids):
        raise ValueError("Mie validation requires all eight shape coefficients to be exactly zero. Regenerate with --sphere.")
    if backend is None:
        try:
            import miepython as backend
        except ImportError as exc:
            raise RuntimeError("Install vis/requirements.txt to enable the Mie reference; no substitute model is used.") from exc
    version = str(getattr(backend, "__version__", "unknown"))
    if not version.startswith("3."):
        raise RuntimeError(f"This adapter requires miepython 3.x (documented API: 3.3.0), got {version}")
    g = view.grid
    a = g.radius_mm
    wavelength_nm = metadata_float(g.metadata, "wavelength_nm", positive=True)
    inside = metadata_float(g.metadata, "interior_index", positive=True)
    outside = metadata_float(g.metadata, "exterior_index", positive=True)
    diameter_nm = 2.0 * a * 1.0e6
    # Supply the absolute sphere index to this high-level API; it handles n_env.
    # Reusing unique theta rows avoids O(Ntheta*Nphi) redundant Mie sums.
    ipar, iper = backend.intensities(complex(inside, 0.0), diameter_nm, wavelength_nm,
                                    np.cos(g.theta), n_env=outside, norm="qsca", n_pole=0)
    ipar = np.asarray(ipar, dtype=np.float64).reshape(-1)
    iper = np.asarray(iper, dtype=np.float64).reshape(-1)
    if ipar.shape != g.theta.shape or iper.shape != g.theta.shape:
        raise RuntimeError("Unexpected miepython return shape")
    if not (np.isfinite(ipar).all() and np.isfinite(iper).all() and (ipar >= 0).all() and (iper >= 0).all()):
        raise RuntimeError("Mie library returned invalid intensities")
    fs, fp = incident_projection_weights(view)
    ds = math.pi * a*a * iper[:, None] * fs
    dp = math.pi * a*a * ipar[:, None] * fp
    result = ds if view.component == "s" else dp if view.component == "p" else ds+dp
    qext, qsca, qback, asymmetry = backend.efficiencies(complex(inside, 0.0), diameter_nm, wavelength_nm, n_env=outside)
    if not np.isfinite([qext, qsca, qback, asymmetry]).all():
        raise RuntimeError("Mie library returned nonfinite efficiencies")
    # This sampled total is an angular-resolution diagnostic, not an exact
    # integral: a sharp forward peak can be unresolved on the supplied grid.
    sampled_total = float(np.sum((ds+dp) * g.solid_angle))
    exact_total = float(math.pi*a*a*qsca)
    info = {
        "miepython_version": version, "backend_JIT": bool(getattr(backend, "USE_JIT", False)),
        "library_normalization": "qsca", "intensities_return_order": ["parallel", "perpendicular"],
        "multipoles": "all (n_pole=0); NOT a Debye/ray-family truncation",
        "radius_mm": a, "diameter_nm": diameter_nm, "vacuum_wavelength_nm": wavelength_nm,
        "absolute_sphere_refractive_index": inside, "environment_refractive_index": outside,
        "relative_refractive_index": inside/outside,
        "size_parameter": 2*math.pi*a*1e6*outside/wavelength_nm,
        "qext": float(qext), "qsca": float(qsca), "qback": float(qback), "g": float(asymmetry),
        "exact_total_cross_section_mm2": exact_total,
        "grid_sampled_total_cross_section_mm2": sampled_total,
        "grid_sampled_to_exact_total_ratio": sampled_total/exact_total if exact_total > 0 else None,
        "sampled_integral_not_guaranteed_angular_convergence": True,
        "comparison_uses_intensity_not_complex_phase": True,
        "mie_angles": "nominal double-precision CSV theta cell centers; C++ query directions are stored in FP32",
    }
    return result, info
