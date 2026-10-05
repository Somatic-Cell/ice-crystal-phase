"""Plotting and explicitly solid-angle-weighted diagnostics; no hidden smoothing."""
from __future__ import annotations

from pathlib import Path
import platform
from typing import Any

import numpy as np
import matplotlib
matplotlib.use("Agg")                    # batch use, including Windows without an IDE
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, BoundaryNorm

from optics_io import OpticsGrid, IntensityView, region_mask


def versions() -> dict[str, str]:
    return {"python": platform.python_version(), "numpy": np.__version__, "matplotlib": matplotlib.__version__}


def plot_heatmap(grid: OpticsGrid, values: np.ndarray, mask: np.ndarray, path: Path,
                 title: str, *, theta_range=(0.0, 180.0), log: bool = True,
                 vmin: float | None = None, vmax: float | None = None,
                 label: str = "Model angular density [mm$^2$/sr]", footer: str = "") -> dict[str, Any]:
    """pcolormesh colors each original cell. No resampling, blur, or log-data export."""
    use = mask & region_mask(grid, theta_range) & np.isfinite(values)
    if log:
        use &= values > 0
    selected = values[use]
    if vmin is not None and (not np.isfinite(vmin) or (log and vmin <= 0)):
        raise ValueError("Invalid vmin")
    if vmax is not None and (not np.isfinite(vmax) or (log and vmax <= 0)):
        raise ValueError("Invalid vmax")
    fig, ax = plt.subplots(figsize=(11.5, 5.8), layout="constrained")
    if selected.size:
        low = float(selected.min()) if vmin is None else vmin
        high = float(selected.max()) if vmax is None else vmax
        if high == low and vmin is None and vmax is None:
            if log: low, high = low / 2, high * 2
            else: low, high = low - 0.5, high + 0.5
        if not low < high:
            raise ValueError("Color limits must satisfy vmin < vmax")
        norm = LogNorm(low, high, clip=False) if log else Normalize(low, high, clip=False)
        masked = np.ma.array(values, mask=~use)
        mesh = ax.pcolormesh(np.linspace(-180, 180, grid.shape[1]+1),
                             np.linspace(0, 180, grid.shape[0]+1), masked,
                             shading="flat", norm=norm, rasterized=True)
        fig.colorbar(mesh, ax=ax, label=label, extend="both")
    else:
        low, high = None, None
        ax.text(0.5, 0.5, "No positive accepted samples" if log else "No accepted samples",
                transform=ax.transAxes, ha="center", va="center")
    ax.set(xlabel="Azimuth phi [deg]", ylabel="Scattering angle theta [deg]",
           xlim=(-180, 180), ylim=(theta_range[1], theta_range[0]), title=title)
    if footer:
        fig.supxlabel(footer, fontsize=8)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return {"log_color_mapping": log, "vmin": low, "vmax": high,
            "positive_or_linear_displayed_samples": int(use.sum()),
            "sample_min": None if not selected.size else float(selected.min()),
            "sample_max": None if not selected.size else float(selected.max()),
            "zero_is_not_replaced_by_epsilon": True}


def plot_status(view: IntensityView, path: Path, *, theta_range=(0.0, 180.0)) -> None:
    # A separate categorical map distinguishes true zero, pending, and numerical
    # failures, even though LogNorm cannot display zero in the intensity map.
    state = np.zeros(view.grid.shape, dtype=np.uint8)
    state[(view.values == 0) & view.known_complete] = 1
    state[~view.known_complete] = 2
    state[~view.numerical_valid | ~np.isfinite(view.values) | (view.values < 0)] = 3
    fig, ax = plt.subplots(figsize=(11.5, 5.8), layout="constrained")
    mesh = ax.pcolormesh(np.linspace(-180, 180, state.shape[1]+1),
                         np.linspace(0, 180, state.shape[0]+1), state,
                         shading="flat", norm=BoundaryNorm(np.arange(-0.5, 4.5), 256), rasterized=True)
    cbar = fig.colorbar(mesh, ax=ax, ticks=[0, 1, 2, 3])
    cbar.ax.set_yticklabels(["known hits complete: positive", "known hits complete: zero",
                            "known hits incomplete", "numerical error / invalid"])
    ax.set(xlabel="Azimuth phi [deg]", ylabel="Scattering angle theta [deg]",
           xlim=(-180, 180), ylim=(theta_range[1], theta_range[0]), title="Evaluation status (NOT full optical validity)")
    fig.supxlabel("Source coverage and omitted physics are global conditions, not encoded by this map.", fontsize=8)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180); plt.close(fig)


def diagnostic_summary(view: IntensityView, theta_range=(0.0, 180.0)) -> dict[str, Any]:
    domain = region_mask(view.grid, theta_range)
    use = view.accepted() & domain
    weight = view.grid.solid_angle
    return {"quantity": "regular_partial_model_density_in_cross_section_units",
            "units": "mm^2/sr", "stage": view.stage, "component": view.component,
            "unpolarized": view.is_unpolarized,
            "theta_range_requested_deg": list(theta_range),
            "selection": "cell centers inside requested range; full cell solid angle is used",
            "grid_shape_theta_phi": list(view.grid.shape),
            "requested_samples": int(domain.sum()), "accepted_known_complete_samples": int(use.sum()),
            "pending_samples": int((domain & view.numerical_valid & ~view.known_complete).sum()),
            "error_samples": int((domain & ~view.numerical_valid).sum()),
            "accepted_zero_samples": int((use & (view.values == 0)).sum()),
            "requested_solid_angle_sr": float(weight[domain].sum()),
            "accepted_solid_angle_sr": float(weight[use].sum()),
            "accepted_solid_angle_fraction": float(weight[use].sum()/weight[domain].sum()),
            "sampled_partial_integral_mm2": float(np.sum(weight[use]*view.values[use])),
            "integral_is_not_certified_total_cross_section": True,
            "warnings": view.warnings(), "versions": versions(),
            "sources": [g.provenance() for g in view.grids]}


def weighted_metrics(model: np.ndarray, reference: np.ndarray, weights: np.ndarray,
                     accepted: np.ndarray) -> dict[str, Any]:
    use = accepted & np.isfinite(model) & np.isfinite(reference) & (model >= 0) & (reference >= 0)
    if not use.any():
        raise ValueError("No common finite accepted samples for comparison")
    a, b, w = model[use], reference[use], weights[use]
    delta = a-b
    denom_l1 = float(np.sum(w*b))
    denom_l2 = float(np.sum(w*b*b))
    positive = (a > 0) & (b > 0)
    lr = np.log10(a[positive]) - np.log10(b[positive])
    return {
        "sample_count": int(use.sum()), "solid_angle_sr": float(w.sum()),
        "model_sampled_integral_mm2": float(np.sum(w*a)),
        "mie_sampled_integral_mm2_same_mask": denom_l1,
        "integral_ratio": float(np.sum(w*a)/denom_l1) if denom_l1 else None,
        "relative_L1_solid_angle_weighted": float(np.sum(w*np.abs(delta))/denom_l1) if denom_l1 else None,
        "relative_L2_solid_angle_weighted": float(np.sqrt(np.sum(w*delta*delta)/denom_l2)) if denom_l2 else None,
        "log10_ratio_rmse_positive_pairs_only": float(np.sqrt(np.sum(w[positive]*lr*lr)/np.sum(w[positive]))) if positive.any() else None,
        "log10_metric_sample_count": int(positive.sum()),
        "zero_model_positive_mie_samples": int(((a == 0) & (b > 0)).sum()),
        "positive_model_zero_mie_samples": int(((a > 0) & (b == 0)).sum()),
        "both_zero_samples": int(((a == 0) & (b == 0)).sum()),
        "scale_fitted": False, "each_curve_independently_normalized": False,
        "log_floor_or_epsilon_added": False,
        "automatic_physical_validation_pass_fail": False,
    }


def plot_slice(grid: OpticsGrid, model: np.ndarray, reference: np.ndarray,
               mask: np.ndarray, path: Path, title: str, phi_deg: float,
               theta_range=(0.0, 180.0), *, mean_azimuth: bool = False) -> dict[str, Any]:
    if mean_azimuth:
        n = mask.sum(axis=1)
        a = np.divide(np.where(mask, model, 0).sum(axis=1), n, out=np.full(n.shape, np.nan), where=n > 0)
        b = np.divide(np.where(mask, reference, 0).sum(axis=1), n, out=np.full(n.shape, np.nan), where=n > 0)
        detail = {"mode": "mean over identical accepted azimuth samples", "accepted_azimuth_samples_per_theta": n.tolist()}
        suffix = "mean over accepted phi samples"
    else:
        delta = (np.degrees(grid.phi)-phi_deg+180) % 360-180
        col = int(np.argmin(np.abs(delta)))
        a = np.where(mask[:, col], model[:, col], np.nan)
        b = np.where(mask[:, col], reference[:, col], np.nan)
        detail = {"requested_phi_deg": phi_deg, "actual_phi_deg": float(np.degrees(grid.phi[col])), "phi_index": col}
        suffix = f"phi={detail['actual_phi_deg']:.6g} deg (no interpolation)"
    # Plot genuine zeros separately in metadata; do not raise them to a log floor.
    a = np.where(a > 0, a, np.nan); b = np.where(b > 0, b, np.nan)
    fig, ax = plt.subplots(figsize=(10.5, 5.3), layout="constrained")
    ax.plot(np.degrees(grid.theta), a, label="C++ regular partial model")
    ax.plot(np.degrees(grid.theta), b, label="Full Lorenz-Mie", linestyle="--")
    if np.isfinite(a).any() or np.isfinite(b).any(): ax.set_yscale("log")
    ax.set(xlim=theta_range, xlabel="Scattering angle theta [deg]", ylabel="Angular density [mm$^2$/sr]",
           title=title+"\n"+suffix)
    ax.legend()
    fig.supxlabel("No focal phase or diffraction is added by Python. Gaps/zeros are not filled.", fontsize=8)
    fig.savefig(path, dpi=180); plt.close(fig)
    return detail
