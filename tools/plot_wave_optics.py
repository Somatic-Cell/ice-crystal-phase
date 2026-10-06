"""Display C++-computed wave stages. No phase correction, filtering, or fitting here.

Separate tool to avoid overwriting the user's locally edited vis/ directory.
Requires NumPy and Matplotlib (already in the previous vis/ requirements).
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import numpy as np


def load(path: Path):
    meta: dict[str, str] = {}
    with path.open(encoding="utf-8-sig") as file:
        for line in file:
            if line.startswith("#"):
                key, sep, value = line[1:].strip().partition("=")
                if sep:
                    if key in meta:
                        raise ValueError(f"Duplicate metadata key: {key}")
                    meta[key] = value
            else:
                break
    if meta.get("format") not in ("rainbow_wave_optics_v1", "rainbow_wave_optics_v2"):
        raise ValueError("Expected rainbow_wave_optics_v1/v2, not the path-only optics CSV.")
    rows, cols = int(meta["theta_count"]), int(meta["phi_count"])
    if rows < 2 or cols < 1:
        raise ValueError("Invalid angular grid.")
    with path.open(encoding="utf-8-sig") as file:
        for header_index, line in enumerate(file):
            if line.strip() and not line.startswith("#"):
                break
        else:
            raise ValueError("CSV header is missing.")
    data = np.genfromtxt(path, delimiter=",", comments="#", names=True,
                         skip_header=header_index, encoding="utf-8-sig")
    return meta, rows, cols, data


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--stage", choices=("incoherent", "path", "focal", "diffraction"), default="diffraction")
    parser.add_argument("--component", choices=("s", "p", "total"), default="total")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--linear", action="store_true")
    parser.add_argument("--vmin", type=float)
    parser.add_argument("--vmax", type=float)
    args = parser.parse_args()
    try:
        meta, rows, cols, data = load(args.csv)
        data = np.atleast_1d(data)
        required = ("direction_id", "theta_index", "phi_index", "theta_rad", "phi_rad",
                    "known_hits_complete", "optical_flags", "focal_valid", "diffraction_valid",
                    f"{args.stage}_{args.component}")
        if any(key not in (data.dtype.names or ()) for key in required):
            raise ValueError("Missing CSV columns.")
        if len(data) != rows * cols:
            raise ValueError("Incomplete grid; refusing to fill missing rows.")
        ids = data["direction_id"]
        if not np.isfinite(ids).all() or not np.equal(ids, np.floor(ids)).all():
            raise ValueError("Invalid direction IDs.")
        order = np.argsort(ids)
        data = data[order]
        if not np.array_equal(data["direction_id"], np.arange(rows * cols)):
            raise ValueError("Duplicate or missing direction IDs.")
        rr, cc = np.divmod(np.arange(rows * cols), cols)
        if not (np.array_equal(data["theta_index"], rr) and np.array_equal(data["phi_index"], cc)):
            raise ValueError("Grid index mismatch.")
        if not (np.allclose(data["theta_rad"], np.pi * (rr + .5) / rows, rtol=0, atol=1e-12)
                and np.allclose(data["phi_rad"], -np.pi + 2*np.pi*(cc+.5)/cols, rtol=0, atol=1e-12)):
            raise ValueError("Angular convention mismatch.")
        radius = float(meta["radius_mm"])
        if meta.get("format") == "rainbow_wave_optics_v2":
            if (meta.get("incident_polarization") != "unpolarized"
                or meta.get("input_states") != "unpolarized_two_orthogonal_unit_Jones_inputs"
                or meta.get("input_coherency") != "0.5,0,0,0.5"):
                raise ValueError("Invalid unpolarized v2 metadata")
            norm2 = float(meta["incident_total_intensity"])
            if norm2 != 1.0:
                raise ValueError("Expected unit incident total intensity")
            state_label = "unpolarized input"
        else:
            jones = np.array([float(x) for x in meta["incident_field"].split(",")])
            norm2 = float(jones @ jones)
            if jones.shape != (4,) or not np.isfinite(jones).all():
                raise ValueError("Invalid input Jones field")
            state_label = "single Jones input"
        if not np.isfinite(norm2) or not norm2 > 0 or not np.isfinite(radius) or radius <= 0:
            raise ValueError("Invalid input intensity or radius.")
        values = data[f"{args.stage}_{args.component}"] * (radius * radius / norm2)
        if args.stage == "focal":
            valid = data["focal_valid"] == 1
        elif args.stage == "diffraction":
            valid = data["diffraction_valid"] == 1
        else:
            valid = (data["known_hits_complete"] == 1) & (data["optical_flags"] == 0)
        accepted = valid & np.isfinite(values) & (values >= 0)
        if not args.linear:
            accepted &= values > 0
        if not accepted.any():
            raise ValueError("No drawable values for this stage; inspect validity flags.")
        vmin = float(np.min(values[accepted])) if args.vmin is None else args.vmin
        vmax = float(np.max(values[accepted])) if args.vmax is None else args.vmax
        if vmin == vmax and args.vmin is None and args.vmax is None:
            vmax = vmin + (abs(vmin) or 1.0)*1e-6
        if not (np.isfinite(vmin) and np.isfinite(vmax) and vmin < vmax and (args.linear or vmin > 0)):
            raise ValueError("Invalid color range.")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm, Normalize
        norm = Normalize(vmin, vmax) if args.linear else LogNorm(vmin, vmax)
        fig, ax = plt.subplots(figsize=(12, 6), layout="constrained")
        image = ax.pcolormesh(np.linspace(-180, 180, cols+1), np.linspace(0, 180, rows+1),
                             np.ma.array(values.reshape(rows, cols), mask=~accepted.reshape(rows, cols)),
                             norm=norm, shading="flat", rasterized=True)
        ax.set_ylim(180, 0)
        ax.set_xlabel("Azimuth phi [deg]")
        ax.set_ylabel("Scattering angle theta [deg]")
        ax.set_title(f"C++ {args.stage} / {args.component}; {state_label}; a={radius:g} mm\n"
                     "Explicit project conventions; source coverage uncertified")
        fig.colorbar(image, ax=ax, label="Model angular density [mm$^2$/sr]")
        args.out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.out, dpi=160)
        plt.close(fig)
        print(f"Saved {args.out}; accepted {int(accepted.sum())}/{len(accepted)}. No holes filled.")
        return 0
    except (OSError, ValueError, KeyError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
