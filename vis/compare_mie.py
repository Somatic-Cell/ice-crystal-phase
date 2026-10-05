"""Quantitative *diagnostic* comparison to full Mie theory on the same grid.

Incomplete model inputs require explicit --allow-partial-model.  This does not
turn pending/error pixels into valid data, nor assert a validation pass.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

import numpy as np
from optics_io import load_optics, make_view, write_json, region_mask, metadata_bool, STAGES, COMPONENTS
from mie_reference import calculate_mie
from analysis_tools import (plot_heatmap, plot_status, plot_slice, diagnostic_summary, weighted_metrics)


def run_comparison(view, out: Path, *, theta_range=(0.0, 180.0), phi_deg=0.0,
                   vmin=None, vmax=None, backend=None) -> dict:
    reference, mie_info = calculate_mie(view, backend=backend)
    g = view.grid
    domain = region_mask(g, theta_range)
    mask = view.accepted() & domain & np.isfinite(reference) & (reference >= 0)
    report = diagnostic_summary(view, theta_range)
    report["mie"] = mie_info
    report["metrics"] = weighted_metrics(view.values, reference, g.solid_angle, mask)
    report["not_a_validation_certificate"] = True
    report["missing_physics_is_not_replaced_in_python"] = True
    report["complete_input_files"] = all(metadata_bool(x.metadata, "optical_complete") and
                                          metadata_bool(x.metadata, "source_coverage_certified") for x in view.grids)
    if not report["complete_input_files"]:
        report["warnings"].append("These metrics include expected model differences; do not interpret them as pure numerical implementation error.")
    if view.stage == "incoherent":
        report["warnings"].append("Non-interfering ray intensity is being compared to full Mie, which includes interference.")
    report["warnings"].append("C++ traces R/TT/TRT/TRRT; the Mie reference includes the full series, including forward diffraction and higher orders.")
    report["warnings"].append("Input values are point samples, not pixel-averaged intensities. Cell-weighted metrics require a resolution convergence study.")
    ratio = mie_info["grid_sampled_to_exact_total_ratio"]
    if ratio is not None and abs(ratio-1) > 0.01:
        report["warnings"].append(f"Mie grid quadrature differs from its exact total by {(ratio-1)*100:.3g}%; forward peaks may be unresolved. No global renormalization was applied.")
    # Common limits for both heatmaps: no independent contrast normalization.
    positives = np.concatenate([view.values[mask & (view.values > 0)], reference[mask & (reference > 0)]])
    if positives.size:
        if vmin is None: vmin = float(positives.min())
        if vmax is None: vmax = float(positives.max())
        if vmin == vmax: vmin, vmax = vmin/2, vmax*2
    out.mkdir(parents=True, exist_ok=True)
    common_footer = "Same units, same accepted pixels, same color limits. Missing physics is not fitted away."
    report["model_display"] = plot_heatmap(g, view.values, mask, out/"model.png", view.label(),
                                            theta_range=theta_range, vmin=vmin, vmax=vmax, footer=common_footer)
    report["mie_display"] = plot_heatmap(g, reference, mask, out/"mie.png", "Full Lorenz-Mie: matching incident polarization",
                                          theta_range=theta_range, vmin=vmin, vmax=vmax, footer=common_footer,
                                          label="Mie differential cross section [mm$^2$/sr]")
    log_ratio = np.full(g.shape, np.nan)
    positive = mask & (view.values > 0) & (reference > 0)
    log_ratio[positive] = np.log10(view.values[positive])-np.log10(reference[positive])
    lim = float(np.max(np.abs(log_ratio[positive]))) if positive.any() else 1.0
    if lim == 0: lim = 1.0
    plot_heatmap(g, log_ratio, positive, out/"log10_ratio.png", "log10(C++ model / full Mie): positive accepted pairs only",
                 theta_range=theta_range, log=False, vmin=-lim, vmax=lim,
                 label="log10 ratio (0 = equal)", footer="No epsilon; zero/reference-zero counts are reported in metrics.json.")
    plot_status(view, out/"status.png", theta_range=theta_range)
    report["slice"] = plot_slice(g, view.values, reference, mask, out/"theta_slice.png", view.label(), phi_deg, theta_range)
    if view.is_unpolarized:
        report["azimuth_mean"] = plot_slice(g, view.values, reference, mask, out/"theta_azimuth_mean.png",
                                              view.label(), phi_deg, theta_range, mean_azimuth=True)
        # Symmetry check only for a true unpolarized pair of sphere calculations.
        full_rows = mask.all(axis=1)
        mean = np.mean(view.values, axis=1)
        good_rows = full_rows & (mean > 0)
        contrast = np.std(view.values[good_rows], axis=1)/mean[good_rows]
        report["unpolarized_axisymmetry"] = {
            "evaluated_fully_accepted_theta_rows": int(good_rows.sum()),
            "phi_std_over_mean_max": float(contrast.max()) if contrast.size else None,
            "phi_std_over_mean_median": float(np.median(contrast)) if contrast.size else None,
            "expect_exact_spherical_unpolarized_intensity_to_be_phi_independent": True,
        }
    else:
        report["warnings"].append("A fixed linear incident polarization generally produces phi-dependent sphere intensity. Do not use axisymmetry of a single x/y run as a validity test.")
    # Full rectangular buffers; invalid/pending values and acceptance remain explicit.
    np.savez_compressed(out/"comparison.npz", theta_rad=g.theta, phi_rad=g.phi, solid_angle_sr=g.solid_angle,
                        model_mm2_per_sr=view.values, mie_mm2_per_sr=reference, accepted=mask,
                        known_hits_complete=view.known_complete, numerical_valid=view.numerical_valid,
                        log10_ratio=log_ratio)
    # Table is restricted to the selected angle range, never silently to valid rows.
    ti, pi = np.nonzero(domain)
    table = np.column_stack((ti, pi, np.degrees(g.theta[ti]), np.degrees(g.phi[pi]),
                             mask[ti, pi], view.values[ti, pi], reference[ti, pi], log_ratio[ti, pi]))
    np.savetxt(out/"comparison.csv", table, delimiter=",", comments="",
               header="theta_index,phi_index,theta_deg,phi_deg,accepted,model_mm2_per_sr,mie_mm2_per_sr,log10_ratio",
               fmt=["%d", "%d", "%.17g", "%.17g", "%d", "%.17g", "%.17g", "%.17g"])
    write_json(out/"metrics.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--orthogonal", type=Path, help="second, orthogonal incident state for an incoherent unpolarized average")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--stage", choices=STAGES, default="path")
    parser.add_argument("--component", choices=COMPONENTS, default="total")
    parser.add_argument("--theta-range", nargs=2, type=float, default=(0.0, 180.0), metavar=("MIN", "MAX"))
    parser.add_argument("--phi-deg", type=float, default=0.0, help="1D slice uses nearest actual phi sample; no resampling")
    parser.add_argument("--vmin", type=float)
    parser.add_argument("--vmax", type=float)
    parser.add_argument("--allow-partial-model", action="store_true",
                        help="acknowledge current model/coverage omissions; output remains diagnostic, not a validation certificate")
    args = parser.parse_args()
    try:
        first = load_optics(args.csv)
        second = load_optics(args.orthogonal) if args.orthogonal else None
        view = make_view(first, args.stage, args.component, second)
        complete = all(metadata_bool(x.metadata, "optical_complete") and
                       metadata_bool(x.metadata, "source_coverage_certified") for x in view.grids)
        if not complete and not args.allow_partial_model:
            raise ValueError("C++ output is an incomplete model. Use --allow-partial-model only for an explicitly diagnostic comparison.")
        if not np.isfinite(args.phi_deg): raise ValueError("phi must be finite")
        result = run_comparison(view, args.out, theta_range=args.theta_range, phi_deg=args.phi_deg,
                                vmin=args.vmin, vmax=args.vmax)
        for warning in result["warnings"]: print("Warning:", warning, file=sys.stderr)
        print(f"Saved diagnostic comparison to {args.out}")
        print("Relative weighted L1:", result["metrics"]["relative_L1_solid_angle_weighted"])
        print("Relative weighted L2:", result["metrics"]["relative_L2_solid_angle_weighted"])
        print("No physical-validation PASS/FAIL was assigned.")
        return 0
    except (ValueError, RuntimeError, OSError, MemoryError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
