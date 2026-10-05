"""Plot existing C++ optical CSV data as an equirectangular, logarithmic map."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from optics_io import load_optics, make_view, write_json, region_mask, STAGES, COMPONENTS
from analysis_tools import plot_heatmap, plot_status, diagnostic_summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="*_optics.csv; not vertices/query/hits CSV")
    parser.add_argument("--orthogonal", type=Path, help="second orthogonal incident Jones state; average intensities")
    parser.add_argument("--stage", choices=STAGES, default="path")
    parser.add_argument("--component", choices=COMPONENTS, default="total")
    parser.add_argument("--out", type=Path, required=True, help="output directory")
    parser.add_argument("--theta-range", nargs=2, type=float, default=(0.0, 180.0), metavar=("MIN", "MAX"))
    parser.add_argument("--linear", action="store_true", help="linear color normalization instead of LogNorm")
    parser.add_argument("--vmin", type=float, help="linear physical density value, even with logarithmic colors")
    parser.add_argument("--vmax", type=float)
    parser.add_argument("--show-partial", action="store_true", help="also DISPLAY finite incomplete-hit partial sums; not for quantitative validation")
    args = parser.parse_args()
    try:
        g = load_optics(args.csv)
        other = load_optics(args.orthogonal) if args.orthogonal else None
        view = make_view(g, args.stage, args.component, other)
        args.out.mkdir(parents=True, exist_ok=True)
        for warning in view.warnings(): print("Warning:", warning, file=sys.stderr)
        use = view.accepted(args.show_partial)
        title = view.label()+f"\na={g.radius_mm:.6g} mm; wavelength={float(g.metadata['wavelength_nm']):.6g} nm"
        footer = "Known-hit mask is not source-coverage certification. No focal-line phase/diffraction is added here."
        display = plot_heatmap(g, view.values, use, args.out/"intensity.png", title,
                               theta_range=args.theta_range, log=not args.linear,
                               vmin=args.vmin, vmax=args.vmax, footer=footer)
        plot_status(view, args.out/"status.png", theta_range=args.theta_range)
        report = diagnostic_summary(view, args.theta_range)
        report["display"] = display
        report["show_incomplete_partial_sums"] = args.show_partial
        report["physical_validation_pass"] = None
        write_json(args.out/"summary.json", report)
        # This cache keeps LINEAR physical data and masks. It is not a colored image
        # nor an independently normalized probability density.
        np.savez_compressed(args.out/"angular_data.npz", theta_rad=g.theta, phi_rad=g.phi,
                            solid_angle_sr=g.solid_angle, density_mm2_per_sr=view.values,
                            known_hits_complete=view.known_complete, numerical_valid=view.numerical_valid,
                            selected_region=region_mask(g, args.theta_range),
                            source_metadata_json=np.array(json.dumps([x.provenance() for x in view.grids])))
        print(f"Saved {args.out / 'intensity.png'}")
        print(f"Saved {args.out / 'status.png'}")
        print(f"Saved {args.out / 'summary.json'} and angular_data.npz")
        return 0
    except (ValueError, RuntimeError, OSError, MemoryError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
