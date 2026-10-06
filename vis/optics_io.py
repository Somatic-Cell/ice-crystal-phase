"""Read the *existing* rainbow_patch_optics_v1 CSV without changing its physics.

Angles are polar scattering angle theta and azimuth phi, NOT Mercator y.
All intensive quantities returned by density() are divided by the incident
Jones norm and multiplied by radius_mm**2: model d(sigma)/dOmega [mm^2/sr].
This unit conversion does not make an incomplete model an exact cross section.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import math
from typing import Any

import numpy as np

FORMAT = "rainbow_patch_optics_v1"
ERROR_MASK = 0x1F
PENDING_MASK = 0xE0
STAGES = ("incoherent", "path")
COMPONENTS = ("s", "p", "total")
INTENSITY_COLUMNS = tuple(f"regular_partial_{s}_{c}" for s in STAGES for c in COMPONENTS)
SELECT_COLUMNS = (
    "direction_id", "theta_index", "phi_index", "theta_rad", "phi_rad", "solid_angle_sr",
    "wx", "wy", "wz", "known_hits_complete", "hit_count", "evaluated_hits", "rejected_hits",
    "refinement_hits", "boundary_hits", "singular_hits", "flags", "query_flags",
) + INTENSITY_COLUMNS


def _require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def metadata_float(meta: dict[str, str], key: str, *, positive: bool = False) -> float:
    try:
        value = float(meta[key])
    except (KeyError, ValueError) as exc:
        raise ValueError(f"Missing/invalid metadata: {key}") from exc
    _require(math.isfinite(value) and (not positive or value > 0), f"Invalid metadata: {key}")
    return value


def metadata_vector(meta: dict[str, str], key: str, count: int) -> np.ndarray:
    try:
        value = np.array([float(x) for x in meta[key].split(",")], dtype=np.float64)
    except (KeyError, ValueError) as exc:
        raise ValueError(f"Missing/invalid metadata: {key}") from exc
    _require(value.shape == (count,) and np.isfinite(value).all(), f"Invalid metadata vector: {key}")
    return value


def metadata_bool(meta: dict[str, str], key: str) -> bool:
    _require(meta.get(key) in ("true", "false"), f"Missing/invalid boolean metadata: {key}")
    return meta[key] == "true"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass
class OpticsGrid:
    path: Path
    metadata: dict[str, str]
    theta: np.ndarray                       # (Ntheta,), radians, increasing
    phi: np.ndarray                         # (Nphi,), radians, increasing
    solid_angle: np.ndarray                 # (Ntheta, Nphi), exact cell measure
    fields: dict[str, np.ndarray]           # matrices, including diagnostics
    sha256: str

    @property
    def shape(self) -> tuple[int, int]:
        return self.solid_angle.shape

    @property
    def radius_mm(self) -> float:
        return metadata_float(self.metadata, "radius_mm", positive=True)

    @property
    def is_unpolarized(self) -> bool:
        return self.metadata.get("incident_polarization") == "unpolarized"

    @property
    def jones(self) -> np.ndarray:
        _require(not self.is_unpolarized, "Unpolarized input has no single Jones vector")
        x = metadata_vector(self.metadata, "incident_field", 4)
        return np.array([complex(x[0], x[1]), complex(x[2], x[3])])

    @property
    def incident_norm2(self) -> float:
        if self.is_unpolarized:
            return metadata_float(self.metadata, "incident_total_intensity", positive=True)
        j = self.jones
        return float(np.vdot(j, j).real)

    @property
    def is_sphere(self) -> bool:
        # Exact zero, not a convenient tolerance that accepts non-spherical drops.
        return bool(np.all(metadata_vector(self.metadata, "coefficients", 8) == 0.0))

    @property
    def numerical_valid(self) -> np.ndarray:
        return (self.fields["flags"].astype(np.uint32) & ERROR_MASK) == 0

    @property
    def known_complete(self) -> np.ndarray:
        return self.fields["known_hits_complete"] == 1

    def density(self, stage: str, component: str) -> np.ndarray:
        _require(stage in STAGES and component in COMPONENTS, "Invalid stage/component")
        return self.fields[f"regular_partial_{stage}_{component}"] * (self.radius_mm**2 / self.incident_norm2)

    def warnings(self) -> list[str]:
        result = []
        for key, message in (
            ("optical_complete", "C++ declares optical_complete=false."),
            ("source_coverage_certified", "Source coverage is NOT certified, even at known_hits_complete=1."),
            ("focal_line_phase_applied", "Focal-line phase is NOT applied."),
            ("diffraction_applied", "Diffraction is NOT applied."),
        ):
            if not metadata_bool(self.metadata, key):
                result.append(message)
        for key in ("missing_source_cells", "no_outgoing_source_cells"):
            n = metadata_float(self.metadata, key)
            if n:
                result.append(f"{key}={int(n)}; this is not a per-direction coverage mask.")
        return result

    def provenance(self) -> dict[str, Any]:
        return {"path": str(self.path.resolve()), "sha256": self.sha256,
                "metadata": self.metadata, "shape_theta_phi": list(self.shape),
                "incident_norm_squared": self.incident_norm2, "warnings": self.warnings()}


def load_optics(path: str | Path) -> OpticsGrid:
    """Strict, index-based loading. Incomplete grids and unsupported schemas fail.

    np.loadtxt reads only required columns; it does not parse the complex field
    columns because visualization/comparison uses already-computed intensities.
    Large final-resolution CSV files still require substantial RAM.
    """
    path = Path(path)
    meta: dict[str, str] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                pair = line[1:].strip().split("=", 1)
                if len(pair) == 2:
                    key, value = (v.strip() for v in pair)
                    _require(key not in meta, f"Duplicate metadata key: {key}")
                    meta[key] = value
                continue
            names = [x.strip() for x in line.split(",")]
            break
        else:
            raise ValueError("CSV has no column header")
        _require(len(set(names)) == len(names), "Duplicate CSV columns")
        missing = set(SELECT_COLUMNS) - set(names)
        _require(not missing, f"Missing optical CSV columns: {sorted(missing)}")
        data = np.loadtxt(stream, delimiter=",", comments="#", dtype=np.float64,
                          usecols=[names.index(n) for n in SELECT_COLUMNS], ndmin=2)
    _require(meta.get("format") in (FORMAT, "rainbow_patch_optics_v2"), f"Expected {FORMAT}; do not pass vertices/patch/query CSV")
    _require(meta.get("order") == "theta_major_phi_minor", "Unsupported/missing angular storage order")
    _require(meta.get("angle_units") == "radians", "Only explicitly declared radians are supported")
    _require(meta.get("quantity") in ("regular_partial_model_angular_density_NOT_phase_function", "model_partial_angular_density_NOT_phase_function"),
             "Unsupported optical quantity; update the adapter explicitly for a new format")
    _require(meta.get("density_units") == "input_field_squared*drop_unit_squared_per_sr",
             "Unsupported density units")
    _require(meta.get("field_components") == "s_perpendicular,p_outgoing_cross_s", "Unsupported output basis")
    if meta.get("format") == "rainbow_patch_optics_v2":
        _require(meta.get("input_states") == "unpolarized_two_orthogonal_unit_Jones_inputs"
                 and meta.get("incident_polarization") == "unpolarized", "Invalid v2 input-state declaration")
        _require(metadata_float(meta, "incident_total_intensity") == 1.0, "Expected unit incident total intensity")
        _require(np.array_equal(metadata_vector(meta, "input_coherency", 4), [0.5, 0, 0, 0.5]), "Unsupported input coherency")
    else:
        _require(meta.get("input_states") == "one_coherent_Jones_state", "Expected a single coherent Jones input per v1 file")
        _require(meta.get("incident_polarization") != "unpolarized", "Unpolarized data needs schema v2")
    _require(not metadata_bool(meta, "normalized"), "Normalized files need a different unit conversion")
    ntheta = metadata_float(meta, "theta_count", positive=True)
    nphi = metadata_float(meta, "phi_count", positive=True)
    _require(ntheta.is_integer() and nphi.is_integer(), "Angular counts must be integers")
    nt, np_ = int(ntheta), int(nphi)
    _require(data.shape == (nt * np_, len(SELECT_COLUMNS)), "Row count does not match the full angular grid")
    col = {n: data[:, i] for i, n in enumerate(SELECT_COLUMNS)}
    integer_names = set(SELECT_COLUMNS[:3]) | set(SELECT_COLUMNS[9:18])
    for name in integer_names:
        x = col[name]
        _require(bool(np.isfinite(x).all() and (x >= 0).all() and (x <= 0xffffffff).all()
                      and (x == np.floor(x)).all()), f"Invalid integer column: {name}")
    ids = col["direction_id"].astype(np.int64)
    expected_ids = np.arange(nt * np_, dtype=np.int64)
    if not np.array_equal(ids, expected_ids):
        idx = np.argsort(ids)
        _require(np.array_equal(ids[idx], expected_ids), "Duplicate or missing direction_id")
        data = data[idx]
        col = {n: data[:, i] for i, n in enumerate(SELECT_COLUMNS)}
    _require(np.array_equal(col["theta_index"], expected_ids // np_), "theta_index/direction_id mismatch")
    _require(np.array_equal(col["phi_index"], expected_ids % np_), "phi_index/direction_id mismatch")
    theta = np.pi * (np.arange(nt) + 0.5) / nt
    phi = -np.pi + 2 * np.pi * (np.arange(np_) + 0.5) / np_
    th = col["theta_rad"].reshape(nt, np_)
    ph = col["phi_rad"].reshape(nt, np_)
    _require(np.allclose(th, theta[:, None], rtol=0, atol=5e-12), "theta is not the declared cell-centered grid")
    _require(np.allclose(ph, phi[None, :], rtol=0, atol=5e-12), "phi is not the declared cell-centered grid")
    # Stable form of dphi * (cos(theta_lo) - cos(theta_hi)).
    weight = np.broadcast_to((4 * np.pi / np_ * np.sin(theta) * np.sin(np.pi / (2 * nt)))[:, None], (nt, np_))
    solid_angle = col["solid_angle_sr"].reshape(nt, np_)
    _require(np.allclose(solid_angle, weight, rtol=2e-8, atol=2e-15), "Incorrect solid-angle cell weights")
    _require(np.isclose(solid_angle.sum(), 4 * np.pi, rtol=1e-9, atol=1e-12), "Grid weights do not sum to 4*pi")
    w = metadata_vector(meta, "incident_direction", 3)
    e0 = metadata_vector(meta, "incident_basis_x", 3)
    _require(abs(np.dot(w, w)-1) < 1e-5 and abs(np.dot(e0, e0)-1) < 1e-5
             and abs(np.dot(w, e0)) < 1e-5, "Invalid incident frame")
    e1 = np.cross(w, e0); e1 /= np.linalg.norm(e1)
    actual = np.column_stack([col[n] for n in ("wx", "wy", "wz")])
    expected = (np.cos(col["theta_rad"])[:, None] * w
                + np.sin(col["theta_rad"])[:, None] * (np.cos(col["phi_rad"])[:, None] * e0
                                                        + np.sin(col["phi_rad"])[:, None] * e1))
    expected /= np.linalg.norm(expected, axis=1)[:, None]
    _require(np.allclose(actual, expected, rtol=0, atol=4e-6), "Directions disagree with angle/frame metadata")
    # Keep the angular intensities and diagnostic matrices, not all input columns.
    fields = {n: col[n].reshape(nt, np_) for n in SELECT_COLUMNS[9:]}
    flags = fields["flags"].astype(np.uint32)
    _require(not np.any(flags & np.uint32(0xffffff00)), "Unknown optical status bits")
    error = (flags & ERROR_MASK) != 0
    expected_complete = ((flags & (ERROR_MASK | PENDING_MASK)) == 0) & (fields["rejected_hits"] == 0) & (fields["evaluated_hits"] == fields["hit_count"])
    _require(np.array_equal(fields["known_hits_complete"], expected_complete), "Completion bit/count inconsistency")
    _require(np.all(fields["evaluated_hits"] <= fields["hit_count"]), "evaluated_hits exceeds hit_count")
    _require(np.all(fields["rejected_hits"] <= fields["hit_count"]), "rejected_hits exceeds hit_count")
    for n in INTENSITY_COLUMNS:
        x = fields[n]
        _require(np.all(np.isfinite(x[~error])) and np.all(x[~error] >= 0), f"Unflagged invalid/negative intensity: {n}")
        _require(np.all(np.isnan(x[error])), f"Numerical-error rows must preserve NaN: {n}")
    for stage in STAGES:
        _require(np.allclose(fields[f"regular_partial_{stage}_total"],
                             fields[f"regular_partial_{stage}_s"] + fields[f"regular_partial_{stage}_p"],
                             rtol=5e-12, atol=0, equal_nan=True), "s+p does not match total")
    grid = OpticsGrid(path, meta, theta, phi, solid_angle, fields, sha256_file(path))
    _require(math.isfinite(grid.incident_norm2) and grid.incident_norm2 > 0, "Zero/nonfinite incident Jones norm")
    _require(np.isclose(metadata_float(meta, "physical_area_factor_mm2", positive=True), grid.radius_mm**2,
                        rtol=5e-12, atol=0), "Physical area scale must be radius_mm squared")
    for name in ("wavelength_nm", "exterior_index", "interior_index"):
        metadata_float(meta, name, positive=True)
    metadata_vector(meta, "coefficients", 8)
    grid.warnings()                         # also validates completion metadata
    return grid


@dataclass
class IntensityView:
    grids: tuple[OpticsGrid, ...]
    stage: str
    component: str
    values: np.ndarray                     # mm^2/sr, incident norm removed
    known_complete: np.ndarray
    numerical_valid: np.ndarray
    is_unpolarized: bool

    @property
    def grid(self) -> OpticsGrid:
        return self.grids[0]

    def accepted(self, show_partial: bool = False) -> np.ndarray:
        return self.numerical_valid & np.isfinite(self.values) & (self.values >= 0) & (True if show_partial else self.known_complete)

    def warnings(self) -> list[str]:
        return list(dict.fromkeys(message for g in self.grids for message in g.warnings()))

    def label(self) -> str:
        state = ("C++ unpolarized input" if self.grid.is_unpolarized else "unpolarized intensity average") if self.is_unpolarized else "single Jones input"
        return f"{self.stage} / {self.component}; {state}; regular partial model"


def make_view(first: OpticsGrid, stage: str = "path", component: str = "total",
              orthogonal: OpticsGrid | None = None) -> IntensityView:
    value = first.density(stage, component)
    if orthogonal is None:
        return IntensityView((first,), stage, component, value, first.known_complete.copy(), first.numerical_valid.copy(), first.is_unpolarized)
    _require(not first.is_unpolarized and not orthogonal.is_unpolarized, "Do not average an already unpolarized file with another input")
    other = orthogonal
    _require(first.shape == other.shape and np.array_equal(first.theta, other.theta)
             and np.array_equal(first.phi, other.phi), "The two input-state grids must match; no resampling is performed")
    # Fields that define particle, experiment, approximation and interpolation.
    numeric_keys = ("radius_mm", "wavelength_nm", "interior_index", "exterior_index", "physical_area_factor_mm2",
                    "incident_grid_half_extent_drop", "reference_distance_drop", "outgoing_reference_distance_drop")
    for key in numeric_keys:
        _require(first.metadata.get(key) == other.metadata.get(key), f"Input-state metadata mismatch: {key}")
    for key in ("coefficients", "incident_direction", "incident_basis_x", "incident_grid", "interpolation", "transport",
                "optical_complete", "source_coverage_certified", "focal_line_phase_applied", "diffraction_applied",
                "missing_source_cells", "no_outgoing_source_cells", "phasor_convention"):
        _require(first.metadata.get(key) == other.metadata.get(key), f"Input-state metadata mismatch: {key}")
    a = first.jones / math.sqrt(first.incident_norm2)
    b = other.jones / math.sqrt(other.incident_norm2)
    _require(abs(np.vdot(a, b)) < 2e-6, "--orthogonal must be an orthogonal Jones state, not a repeated run")
    # Do not average fields. Normalize each incident state before intensity averaging.
    value = 0.5 * (value + other.density(stage, component))
    return IntensityView((first, other), stage, component, value,
                         first.known_complete & other.known_complete,
                         first.numerical_valid & other.numerical_valid, True)


def region_mask(grid: OpticsGrid, theta_range: tuple[float, float] | list[float]) -> np.ndarray:
    lo, hi = theta_range
    _require(0 <= lo < hi <= 180, "theta range must satisfy 0 <= lo < hi <= 180")
    rows = (np.degrees(grid.theta) >= lo) & (np.degrees(grid.theta) <= hi)
    _require(bool(rows.any()), "No cell centers in the requested theta range")
    return np.broadcast_to(rows[:, None], grid.shape)


def write_json(path: Path, value: dict[str, Any]) -> None:
    def clean(x: Any) -> Any:
        if isinstance(x, dict): return {str(k): clean(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)): return [clean(v) for v in x]
        if isinstance(x, np.ndarray): return clean(x.tolist())
        if isinstance(x, np.generic): return clean(x.item())
        if isinstance(x, float) and not math.isfinite(x): return None
        if isinstance(x, Path): return str(x)
        return x
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean(value), ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
