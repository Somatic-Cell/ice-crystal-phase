#!/usr/bin/env python3
"""Inspect a saved rainbow CDF without tracing, sampling, fitting, or filtering.

Place in <repository>/tools/. Requires the existing python/rainbow_dataset/reader.
NumPy + Matplotlib are required; CUDA-enabled PyTorch is optional. All CDF
arithmetic is binary64. The PNG contains log10(p_omega / (1 sr^-1)).
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import importlib
import json
import math
from pathlib import Path
import sys
import time
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if (ROOT / "python").is_dir():
    sys.path.insert(0, str(ROOT / "python"))


@dataclass
class PlotData:
    # Display order: theta rows, phi columns. This owns its memory, not an mmap.
    log10_pdf: np.ndarray
    theta_edges_deg: np.ndarray
    phi_edges_deg: np.ndarray
    theta_range: tuple[float, float]
    phi_range: tuple[float, float]
    summary: dict[str, Any]
    metadata: dict[str, Any]


class Backend:
    """Use NumPy or PyTorch tensors. No CPU fallback after a CUDA error."""

    def __init__(self, requested: str) -> None:
        self.torch: Any = None
        self.device: Any = None
        self.name = "numpy/cpu"
        self.selection_note = "CPU explicitly selected"
        if requested == "cpu":
            return
        if requested != "auto" and requested != "cuda" and not requested.startswith("cuda:"):
            raise ValueError("--device must be auto, cpu, cuda, or cuda:N")
        try:
            torch = importlib.import_module("torch")
        except ImportError as error:
            if requested == "auto":
                self.selection_note = "PyTorch not installed; auto selected NumPy/CPU"
                return
            raise RuntimeError("CUDA requested but PyTorch cannot be imported") from error
        if not torch.cuda.is_available():
            if requested == "auto":
                self.selection_note = "CUDA unavailable; auto selected NumPy/CPU"
                return
            raise RuntimeError("CUDA requested but torch.cuda.is_available() is False")
        device = torch.device("cuda:0" if requested in ("auto", "cuda") else requested)
        if device.index is None or not 0 <= device.index < torch.cuda.device_count():
            raise ValueError("Invalid CUDA device ordinal")
        self.torch, self.device = torch, device
        self.name = f"torch/{device} ({torch.cuda.get_device_name(device)})"
        self.selection_note = "CUDA selected; reconstruction, validation and reductions use FP64 tensors"

    def array(self, value: np.ndarray) -> Any:
        if self.torch is None:
            return np.asarray(value, dtype=np.float64)
        # Never expose a read-only mmap as a writable torch.from_numpy tensor.
        # This is a bounded staging copy, not a copy of the entire CDF table.
        owned = np.array(value, dtype=np.float64, order="C", copy=True)
        return self.torch.from_numpy(owned).to(device=self.device, dtype=self.torch.float64)

    def check_block(self, c: Any, first_phi: int) -> None:
        if self.torch is None:
            checks = [
                np.isfinite(c).all(), ((c >= 0) & (c <= 1)).all(),
                (c[:, 0] == 0).all(), (c[:, -1] == 1).all(),
                (c[:, 1:] >= c[:, :-1]).all(),
            ]
        else:
            t = self.torch
            # One small readback per block, not one synchronization per element.
            checks = t.stack([
                t.isfinite(c).all(), ((c >= 0) & (c <= 1)).all(),
                (c[:, 0] == 0).all(), (c[:, -1] == 1).all(),
                (c[:, 1:] >= c[:, :-1]).all(),
            ]).cpu().tolist()
        labels = ("nonfinite data", "out-of-range data", "nonzero first endpoint",
                  "last endpoint not one", "decreasing CDF")
        for okay, label in zip(checks, labels):
            if not okay:
                raise ValueError(f"Conditional CDF: {label} in phi rows starting at {first_phi}; no repair performed")

    def evaluate(self, c: Any, q: Any, area: Any, mean_mu: Any) -> tuple[Any, list[float]]:
        dc = c[:, 1:] - c[:, :-1]
        mass = dc * q[:, None]
        # Factorized logarithm prevents a positive q*dc from turning into an
        # artificial zero solely because its FP64 product underflows.
        if self.torch is None:
            with np.errstate(divide="ignore", invalid="raise"):
                logp = np.log10(dc) + np.log10(q[:, None]) - np.log10(area[None, :])
            zero = (dc == 0) | (q[:, None] == 0)
            underflow = (~zero) & (mass == 0)
            positive = np.isfinite(logp)
            summary = [float(mass.sum(dtype=np.float64)),
                       float((mass * mean_mu[None, :]).sum(dtype=np.float64)),
                       float(zero.sum()), float(underflow.sum()),
                       float(np.min(logp, where=positive, initial=np.inf)),
                       float(np.max(logp, where=positive, initial=-np.inf))]
            if np.isnan(logp).any() or np.isposinf(logp).any():
                raise ValueError("Nonfinite reconstructed log-PDF")
        else:
            t = self.torch
            logp = t.log10(dc) + t.log10(q[:, None]) - t.log10(area[None, :])
            zero = (dc == 0) | (q[:, None] == 0)
            underflow = (~zero) & (mass == 0)
            positive = t.isfinite(logp)
            values = t.stack([
                mass.sum(), (mass * mean_mu[None, :]).sum(),
                zero.sum(dtype=t.float64), underflow.sum(dtype=t.float64),
                t.where(positive, logp, math.inf).amin(),
                t.where(positive, logp, -math.inf).amax(),
                (t.isnan(logp) | t.isposinf(logp)).sum(dtype=t.float64),
            ]).cpu().tolist()
            if values[-1] != 0:
                raise ValueError("Nonfinite reconstructed log-PDF")
            summary = values[:-1]
        return logp, summary

    def download(self, value: Any) -> np.ndarray:
        return np.asarray(value) if self.torch is None else value.cpu().numpy()


def _range(value: tuple[float, float], lower: float, upper: float, name: str) -> tuple[float, float]:
    a, b = map(float, value)
    if not (math.isfinite(a) and math.isfinite(b) and lower <= a < b <= upper):
        raise ValueError(f"{name} must satisfy {lower} <= LOW < HIGH <= {upper}")
    return a, b


def _overlap(edges: np.ndarray, limits: tuple[float, float]) -> tuple[int, int]:
    # Include cells intersecting the view. Do not resample or interpolate values.
    start = max(0, int(np.searchsorted(edges, limits[0], side="right")) - 1)
    stop = min(len(edges) - 1, int(np.searchsorted(edges, limits[1], side="left")))
    if stop <= start:
        raise ValueError("Requested view contains no cells")
    return start, stop


def open_record(directory: Path) -> Any:
    try:
        from rainbow_dataset import PhaseRecord
    except ImportError as error:
        raise RuntimeError(
            "Cannot import the existing PhaseRecord reader. Place this script in tools/ "
            "beside python/rainbow_dataset/, or set PYTHONPATH to that python/ directory."
        ) from error
    # File layout/frame checks remain in the reader. The large CDF scan and HG
    # check are done once below, using the selected backend; they are NOT skipped.
    return PhaseRecord(directory, validate=False)


def reconstruct(record: Any, *, backend: Backend, chunk_phi: int = 256,
                theta_range: tuple[float, float] = (0.0, 180.0),
                phi_range: tuple[float, float] = (-180.0, 180.0)) -> PlotData:
    if type(chunk_phi) is not int or chunk_phi < 1:
        raise ValueError("chunk_phi must be a positive integer")
    tr = _range(theta_range, 0.0, 180.0, "theta range")
    pr = _range(phi_range, -180.0, 180.0, "phi range")
    m = record.metadata
    if m.get("theta_definition") != "angle_between_incident_and_outgoing_propagation_directions":
        raise ValueError("Unsupported theta convention")
    phi_limits = np.asarray(m.get("phi_range_rad"), dtype=np.float64)
    if phi_limits.shape != (2,) or not np.allclose(phi_limits, [-np.pi, np.pi], rtol=0, atol=1e-14):
        raise ValueError("Unsupported phi range")
    u = np.array(record.u_edges, dtype=np.float64, copy=True)
    marginal = np.array(record.phi_cdf, dtype=np.float64, copy=True)
    if (not np.isfinite(u).all() or u[0] != 0 or u[-1] != 1
            or not (np.diff(u) > 0).all()):
        raise ValueError("Invalid u_edges; no repair performed")
    if (not np.isfinite(marginal).all() or marginal[0] != 0 or marginal[-1] != 1
            or not (np.diff(marginal) >= 0).all()):
        raise ValueError("Invalid phi CDF; no repair performed")
    area = (4.0 * np.pi / record.np) * np.diff(u)
    if not np.isfinite(area).all() or not (area > 0).all():
        raise ValueError("Invalid or unrepresentable cell solid angle")
    # Stable at both poles; use SAVED u_edges, not freshly chosen angular nodes.
    theta_edges = np.degrees(2 * np.arctan2(np.sqrt(u), np.sqrt(1.0 - u)))
    phi_edges = np.linspace(-180.0, 180.0, record.np + 1, dtype=np.float64)
    i0, i1 = _overlap(theta_edges, tr)
    j0, j1 = _overlap(phi_edges, pr)
    # Phi-major contiguous storage while processing; transpose only the view.
    view = np.empty((j1 - j0, i1 - i0), dtype=np.float64)
    q = backend.array(np.diff(marginal))
    area_b = backend.array(area)
    mu_b = backend.array(1.0 - (u[:-1] + u[1:]))
    totals: list[float] = []
    moments: list[float] = []
    zero_count = underflow_count = 0
    global_min, global_max = math.inf, -math.inf
    started = time.perf_counter()
    for first in range(0, record.np, chunk_phi):
        stop = min(first + chunk_phi, record.np)
        c = backend.array(record.theta_cdf[first:stop])
        backend.check_block(c, first)
        logp, values = backend.evaluate(c, q[first:stop], area_b, mu_b)
        total, moment, zeros, underflows, low, high = values
        totals.append(total)
        moments.append(moment)
        zero_count += int(zeros)
        underflow_count += int(underflows)
        global_min, global_max = min(global_min, low), max(global_max, high)
        left, right = max(first, j0), min(stop, j1)
        if left < right:
            view[left - j0:right - j0] = backend.download(logp[left - first:right - first, i0:i1])
        del c, logp
    # Scalar summaries only. CUDA work is synchronized by the above downloads.
    mass, moment = math.fsum(totals), math.fsum(moments)
    if not math.isfinite(mass) or abs(mass - 1.0) > 1e-10:
        raise ValueError(f"Reconstructed mass is {mass:.17g}, not one; refusing to renormalize")
    g = moment / mass
    if not math.isfinite(g) or not -1.0 <= g <= 1.0:
        raise ValueError("Invalid reconstructed first moment")
    metadata_g: float | None = None
    if m["schema"] == "rainbow.phase_cdf.numpy.v2":
        hg = m.get("hg", {})
        if (not isinstance(hg, dict) or hg.get("method") != "first_moment_of_saved_cell_pdf"
                or hg.get("target") != "saved_cdf"
                or hg.get("cosine_convention") != "dot(incident_propagation,outgoing_propagation)"):
            raise ValueError("Unsupported HG metadata convention")
        label = hg.get("g")
        if type(label) not in (float, int) or not math.isfinite(label) or not -1 < label < 1:
            raise ValueError("Invalid HG metadata value")
        metadata_g = float(label)
        if abs(g - metadata_g) > 5e-12:
            raise ValueError(f"HG label mismatch: reconstructed={g:.17g}, metadata={metadata_g:.17g}")
    if not math.isfinite(global_min) or not math.isfinite(global_max):
        raise ValueError("No positive finite PDF values")
    summary = {
        "backend": backend.name, "backend_note": backend.selection_note,
        "saved_grid_theta_phi": [record.nt, record.np],
        "cdf_reconstructed_mass": mass, "g_reconstructed": g, "g_metadata": metadata_g,
        "global_log10_pdf_min_positive": global_min, "global_log10_pdf_max": global_max,
        "zero_probability_cells": zero_count,
        "positive_mass_product_underflows": underflow_count,
        "display_theta_indices_half_open": [i0, i1], "display_phi_indices_half_open": [j0, j1],
        "display_shape_theta_phi": [i1 - i0, j1 - j0],
        "display_zero_cells": int(np.isneginf(view).sum()),
        "theta_range_deg": list(tr), "phi_range_deg": list(pr),
        "chunk_phi": chunk_phi,
        "reconstruction_and_validation_seconds": time.perf_counter() - started,
        "pdf_measure": "solid_angle_sr", "cdf_repaired": False,
        "renormalized": False, "smoothed": False,
        "log_definition": "log10(p_omega / (1 sr^-1)); -inf for exact zero cells",
        "source_quality": m.get("quality", {}),
        "storage_processing": m.get("storage_processing", {}),
    }
    return PlotData(view.T, theta_edges[i0:i1 + 1].copy(), phi_edges[j0:j1 + 1].copy(),
                    tr, pr, summary, dict(m))


def plot_limits(data: PlotData, explicit: tuple[float, float] | None,
                decades: float | None) -> tuple[float, float]:
    if explicit is not None:
        lo, hi = map(float, explicit)
        if not math.isfinite(lo) or not math.isfinite(hi) or not lo < hi:
            raise ValueError("--log-limits requires finite LOW < HIGH")
        return lo, hi
    hi = data.summary["global_log10_pdf_max"]
    lo = data.summary["global_log10_pdf_min_positive"]
    # Color limits always refer to the GLOBAL PDF, even for a crop.
    if decades is not None:
        if not math.isfinite(decades) or decades <= 0:
            raise ValueError("--decades must be positive and finite")
        return hi - decades, hi
    if hi - lo < 1e-8:
        center = 0.5 * (lo + hi)
        return center - 0.5, center + 0.5
    return lo, hi


def _uniform(edges: np.ndarray) -> bool:
    return bool(np.allclose(edges, np.linspace(edges[0], edges[-1], len(edges)),
                            rtol=0, atol=1e-10))


def draw(data: PlotData, output: Path, *, log_limits: tuple[float, float], dpi: int = 180,
         cell_image: Path | None = None, show: bool = False) -> None:
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    if dpi < 1:
        raise ValueError("DPI must be positive")
    lo, hi = log_limits
    # Only true zero cells are masked. Corrupt input never reaches drawing.
    image = np.ma.array(data.log10_pdf, mask=np.isneginf(data.log10_pdf), copy=False)
    norm = Normalize(vmin=lo, vmax=hi, clip=False)
    fig = plt.figure(figsize=(12.8, 7.6), dpi=dpi)
    ax = fig.add_axes((0.075, 0.15, 0.78, 0.70))
    cax = fig.add_axes((0.89, 0.15, 0.018, 0.70))
    try:
        if _uniform(data.theta_edges_deg) and _uniform(data.phi_edges_deg):
            artist = ax.imshow(image, origin="upper", aspect="auto", interpolation="nearest",
                               resample=False, norm=norm,
                               extent=(data.phi_edges_deg[0], data.phi_edges_deg[-1],
                                       data.theta_edges_deg[-1], data.theta_edges_deg[0]))
            renderer = "imshow_nearest_no_resample"
        else:
            artist = ax.pcolormesh(data.phi_edges_deg, data.theta_edges_deg, image,
                                   shading="flat", antialiased=False, rasterized=True, norm=norm)
            renderer = "pcolormesh_flat_no_antialias"
        ax.set_xlim(*data.phi_range)
        ax.set_ylim(data.theta_range[1], data.theta_range[0])
        ax.set_xlabel("Azimuth phi [deg]")
        ax.set_ylabel("Scattering angle theta [deg]")
        m = data.metadata
        detail = [f"stage={m.get('stage', 'unknown')}"]
        for key, label in (("radius_mm", "radius [mm]"), ("wavelength_nm", "wavelength [nm]"),
                           ("incident_inclination_degrees", "inclination [deg]")):
            if isinstance(m.get(key), (int, float)):
                detail.append(f"{label}={m[key]:.6g}")
        ax.set_title("Saved-CDF reconstructed PDF; unpolarized input\n" + "; ".join(detail), pad=13)
        visible = data.log10_pdf[np.isfinite(data.log10_pdf)]
        below = bool(visible.size and visible.min() < lo)
        above = bool(visible.size and visible.max() > hi)
        extend = "both" if below and above else "min" if below else "max" if above else "neither"
        cb = fig.colorbar(artist, cax=cax, extend=extend)
        cb.set_label("log10(p_omega / (1 sr^-1))")
        s = data.summary
        footer = (f"saved grid={s['saved_grid_theta_phi'][0]} x {s['saved_grid_theta_phi'][1]}; "
                  f"mass={s['cdf_reconstructed_mass']:.12g}; g={s['g_reconstructed']:.9g}\n"
                  f"Blank cells: p=0 (view={s['display_zero_cells']}, global={s['zero_probability_cells']}). "
                  "No smoothing, CDF repair, or renormalization.")
        quality = m.get("quality", {})
        if isinstance(quality, dict) and quality.get("underresolved_directions", 0):
            footer += f"\nGenerator warning: {quality['underresolved_directions']} underresolved directions."
        fig.text(0.5, 0.038, footer, ha="center", va="bottom", fontsize=9)
        fig.canvas.draw()
        pixels = ax.get_window_extent().size
        s["render_method"] = renderer
        s["axes_pixels_width_height"] = [float(x) for x in pixels]
        s["color_log10_limits"] = [lo, hi]
        s["color_saturated_below_cells_in_view"] = int(np.count_nonzero(visible < lo))
        s["color_saturated_above_cells_in_view"] = int(np.count_nonzero(visible > hi))
        nt, np_ = data.log10_pdf.shape
        if np_ > pixels[0] or nt > pixels[1]:
            print("[cdf-plot] WARNING: more CDF cells than figure pixels; screen/PNG resampling can hide "
                  "detail. Use a cropped view, larger --dpi, or --cell-image and inspect it at 100%.",
                  file=sys.stderr)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=dpi)
        if cell_image is not None:
            cell_image.parent.mkdir(parents=True, exist_ok=True)
            # One stored cell -> one RGBA pixel, NO resizing. Index-space image:
            # irregular theta cells do not have their angular widths represented.
            plt.imsave(cell_image, image, vmin=lo, vmax=hi, origin="upper")
        if show:
            plt.show()
    finally:
        plt.close(fig)


def _output_path(path: Path, *, record_dir: Path, overwrite: bool, suffix: str) -> Path:
    resolved = path.resolve()
    if resolved.suffix.lower() != suffix:
        raise ValueError(f"Output must have {suffix} suffix: {path}")
    # Keep generated plots separate from the immutable training record.
    if resolved.is_relative_to(record_dir.resolve()):
        raise ValueError("Write debug outputs outside the dataset-record directory")
    if resolved.exists() and not overwrite:
        raise FileExistsError(f"Output exists: {path}; specify --overwrite to replace it")
    return resolved


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("record", type=Path, help="directory containing metadata.json and three NPY arrays")
    parser.add_argument("--out", required=True, type=Path, help="output PNG outside the record directory")
    parser.add_argument("--device", default="auto", help="auto (default), cpu, cuda, cuda:N")
    parser.add_argument("--chunk-phi", type=int, default=256, help="conditional CDF rows per transfer")
    parser.add_argument("--theta-range", nargs=2, type=float, default=(0.0, 180.0), metavar=("LOW", "HIGH"))
    parser.add_argument("--phi-range", nargs=2, type=float, default=(-180.0, 180.0), metavar=("LOW", "HIGH"))
    color = parser.add_mutually_exclusive_group()
    color.add_argument("--log-limits", nargs=2, type=float, metavar=("LOW", "HIGH"), help="fixed log10 color limits")
    color.add_argument("--decades", type=float, help="color range below the GLOBAL maximum; PDF is unchanged")
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--cell-image", type=Path, help="optional one-cell-per-pixel PNG without axes")
    parser.add_argument("--report", type=Path, help="optional JSON summary")
    parser.add_argument("--show", action="store_true", help="also open the Matplotlib window")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    try:
        out = _output_path(args.out, record_dir=args.record, overwrite=args.overwrite, suffix=".png")
        raw = None if args.cell_image is None else _output_path(
            args.cell_image, record_dir=args.record, overwrite=args.overwrite, suffix=".png")
        report = None if args.report is None else _output_path(
            args.report, record_dir=args.record, overwrite=args.overwrite, suffix=".json")
        destinations = [p for p in (out, raw, report) if p is not None]
        if len(set(destinations)) != len(destinations):
            raise ValueError("Output paths must be distinct")
        backend = Backend(args.device)
        print(f"[cdf-plot] {backend.name}: {backend.selection_note}")
        with open_record(args.record) as record:
            data = reconstruct(record, backend=backend, chunk_phi=args.chunk_phi,
                               theta_range=tuple(args.theta_range), phi_range=tuple(args.phi_range))
        limits = plot_limits(data, None if args.log_limits is None else tuple(args.log_limits), args.decades)
        draw(data, out, log_limits=limits, dpi=args.dpi, cell_image=raw, show=args.show)
        if report is not None:
            report.parent.mkdir(parents=True, exist_ok=True)
            with report.open("w", encoding="utf-8") as stream:
                json.dump(data.summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
                stream.write("\n")
        s = data.summary
        print(f"[cdf-plot] mass={s['cdf_reconstructed_mass']:.17g}, g={s['g_reconstructed']:.17g}, "
              f"metadata_g={s['g_metadata']}, zero_cells={s['zero_probability_cells']}, "
              f"positive_mass_product_underflows={s['positive_mass_product_underflows']}")
        print(f"[cdf-plot] log10 PDF range=[{s['global_log10_pdf_min_positive']:.8g}, "
              f"{s['global_log10_pdf_max']:.8g}], color range={limits}; saved {out}")
        return 0
    except (OSError, ValueError, RuntimeError, ImportError, KeyError, TypeError) as error:
        print(f"plot_phase_cdf: {type(error).__name__}: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
