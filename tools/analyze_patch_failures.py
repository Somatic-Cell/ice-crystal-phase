"""Classify recorded patch witnesses without changing simulator outputs.

Uses exact rational arithmetic for J=P dot (P_u cross P_v) at the four
corners of the STORED binary32 bilinear patch. This is not an exact physical
wavefront or raindrop-intersection solver. No epsilon makes a zero positive.
Python >=3.10; standard library only.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import struct
import sys
from typing import Any

FORMAT = "rainbow_patch_failure_witness_v1"
FAMILIES = ("R", "TT", "TRT", "TRRT")


def integer(value: Any, name: str, maximum: int = 0xFFFFFFFF) -> int:
    if type(value) is not int or not 0 <= value <= maximum:
        raise ValueError(f"{name}: invalid unsigned integer")
    return value


def f32(bits: int) -> float:
    return struct.unpack(">f", struct.pack(">I", integer(bits, "float bits")))[0]


def rational(bits: int) -> Fraction:
    value = f32(bits)
    if not math.isfinite(value):
        raise ValueError("Nonfinite input has no exact rational representation")
    # Every binary32 is exactly representable in Python's binary64 float.
    return Fraction.from_float(value)


def vector(bits: list[int]) -> tuple[Fraction, ...]:
    if len(bits) != 3:
        raise ValueError("Expected three vector components")
    return tuple(rational(x) for x in bits)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])


def sign(x: Fraction) -> int:
    return (x > 0) - (x < 0)


def exact_patch(patch: dict) -> dict:
    vertices = patch["vertices"]
    if len(vertices) != 4:
        raise ValueError("Patch must have four vertices in 00,10,01,11 order")
    try:
        a, b, c, d = (vector(v["direction_bits"]) for v in vertices)
    except ValueError as exc:
        return {"classification": "nonfinite_corner", "detail": str(exc), "signs": [], "jacobians": []}
    j = (dot(a, cross(sub(b, a), sub(c, a))),
         dot(b, cross(sub(b, a), sub(d, b))),
         dot(c, cross(sub(d, c), sub(c, a))),
         dot(d, cross(sub(d, c), sub(d, b))))
    signs = [sign(x) for x in j]
    axis_bits = patch["host_audit"]["hemisphere_axis_bits"]
    try:
        axis = vector(axis_bits)
        hemisphere = [dot(axis, x) for x in (a, b, c, d)]
        original_hemisphere = all(x > 0 for x in hemisphere)
    except ValueError:
        hemisphere = []
        original_hemisphere = False
    # Alternative exact-sum witness is diagnostic only, not the production axis.
    alternate_axis = add(add(a, b), add(c, d))
    alternate_h = [dot(alternate_axis, x) for x in (a, b, c, d)]
    alternate_hemisphere = all(x > 0 for x in alternate_h)
    if 1 in signs and -1 in signs:
        kind = "jacobian_sign_change"
    elif 0 in signs:
        kind = "jacobian_zero_at_corner"
    elif not original_hemisphere:
        kind = "original_hemisphere_candidate_not_certified"
    else:
        kind = "regular_positive" if signs[0] > 0 else "regular_negative"
    return {"classification": kind, "signs": signs,
            "jacobians": [str(x) for x in j],
            "jacobian_values": [float(x) for x in j],
            "original_hemisphere_certified": original_hemisphere,
            "hemisphere_dots": [str(x) for x in hemisphere],
            "alternate_exact_sum_hemisphere_certified": alternate_hemisphere,
            "alternate_hemisphere_dots": [str(x) for x in alternate_h]}


def validate_enclosures(exact: dict, patch: dict) -> list[str]:
    """Audit host interval outputs against rationals; never assumes FP64 exact."""
    issues: list[str] = []
    for key, reference in (("jacobian", exact.get("jacobians", [])),
                           ("hemisphere", exact.get("hemisphere_dots", []))):
        if not reference:
            continue
        for width in (32, 64):
            bounds = patch["host_audit"][f"{key}{width}"]
            if len(bounds) != 4:
                raise ValueError("Wrong interval count")
            for k, (pair, value) in enumerate(zip(bounds, reference)):
                if len(pair) != 2:
                    raise ValueError("Wrong interval shape")
                lo, hi = pair
                if lo is None or hi is None:
                    issues.append(f"{key}{width}[{k}]: nonfinite bound")
                elif not Fraction.from_float(float(lo)) <= Fraction(value) <= Fraction.from_float(float(hi)):
                    issues.append(f"{key}{width}[{k}]: exact value outside interval")
    return issues


def classify_reason(patch: dict, exact: dict) -> str:
    kind = exact["classification"]
    if kind == "jacobian_sign_change":
        return "true_sign_change_in_stored_map"
    if kind == "jacobian_zero_at_corner":
        return "exact_zero_not_repairable_by_more_precision"
    if kind == "nonfinite_corner":
        return "invalid_corner_data"
    if kind == "original_hemisphere_candidate_not_certified":
        return "hemisphere_candidate_issue_NOT_proof_no_hemisphere_exists"
    a = patch["host_audit"]
    direction = exact["signs"][0]
    omega = f32(a["signed_omega32_bits"])
    if a["orientation32"] == 0:
        return ("fp32_certificate_inconclusive_fp64_certifies" if a["interval_orientation64"] == direction
                else "intervals_inconclusive_exact_signs_certify")
    if not math.isfinite(omega) or omega * direction <= 0:
        return "signed_solid_angle_metric_issue"
    if patch["stored_status"] == 3:
        return "stored_pending_but_host_predicate_and_metric_pass"
    if patch["stored_status"] not in (1, 2):
        return "unexpected_nonregular_stored_status"
    if patch["stored_status"] != (1 if direction > 0 else 2):
        return "stored_orientation_disagrees_with_exact_sign"
    return "regular_context_patch"


def no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def load(path: Path) -> dict:
    with path.open(encoding="utf-8-sig") as f:
        report = json.load(f, object_pairs_hook=no_duplicates,
                           parse_constant=lambda x: (_ for _ in ()).throw(ValueError(f"Invalid JSON {x}")))
    if report.get("format") != FORMAT or report.get("report_complete") is not True:
        raise ValueError("Wrong format or incomplete .part report")
    if report.get("changes_numerical_results") is not False:
        raise ValueError("Expected a read-only witness report")
    return report


def first_stage(d: dict) -> str:
    # The masks are from baseline 8c8deab. Benign FP64/underresolved bits excluded.
    if d["query_flags"] & (1 | 2 | 4 | 256):
        return "query_error"
    if d["optical_flags"] or d["rejected_hits"] or d["evaluated_hits"] != d["hit_count"]:
        return "optical"
    if d["focal_flags"]:
        return "focal"
    if d["diffraction_flags"] & (4 | 8 | 16 | 32):
        return "diffraction"
    return "nonfinite_without_explanatory_flag"


def csv_write(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def analyze(report: dict) -> tuple[dict, list[dict], list[dict], list[dict]]:
    config = report["config"]
    nx, ny = integer(config["grid_width"], "grid_width") - 1, integer(config["grid_height"], "grid_height") - 1
    nt, np_ = integer(config["theta_count"], "theta_count"), integer(config["phi_count"], "phi_count")
    if min(nx, ny, nt, np_) < 1:
        raise ValueError("Invalid grid")
    cells = nx * ny
    patch_rows, direction_rows, hit_rows = [], [], []
    exact_by_id, patch_by_compact, reasons = {}, {}, {}
    all_issues = []
    for patch in report["patches"]:
        pid = integer(patch["patch_id"], "patch_id", 4*cells-1)
        ci = integer(patch["compact_index"], "compact_index")
        if pid in exact_by_id or ci in patch_by_compact:
            raise ValueError("Duplicate patch identity")
        exact = exact_patch(patch)
        reason = classify_reason(patch, exact)
        issues = validate_enclosures(exact, patch)
        for issue in issues:
            all_issues.append({"patch_id": pid, "issue": issue})
        family, cell = divmod(pid, cells)
        iy, ix = divmod(cell, nx)
        if (family, ix, iy) != (patch["family"], patch["cell_x"], patch["cell_y"]):
            raise ValueError("Inconsistent patch/cell identity")
        expected_first = family*(nx+1)*(ny+1) + iy*(nx+1) + ix
        expected = [expected_first, expected_first+1, expected_first+nx+1, expected_first+nx+2]
        if [v["vertex_index"] for v in patch["vertices"]] != expected:
            all_issues.append({"patch_id": pid, "issue": "vertex indices disagree with uniform-grid layout"})
        row = {"patch_id": pid, "family": FAMILIES[family], "cell_x": ix, "cell_y": iy,
               "stored_status": patch["stored_status"], "classification": exact["classification"],
               "reason": reason, "jacobian_signs": ",".join(map(str, exact["signs"])),
               "host_orientation32": patch["host_audit"]["orientation32"],
               "interval_orientation64": patch["host_audit"]["interval_orientation64"],
               "center_exactly_on_axis": (2*ix+1 == nx and 2*iy+1 == ny),
               "interval_audit_issues": len(issues)}
        patch_rows.append(row)
        exact_by_id[pid] = exact
        patch_by_compact[ci] = patch
        reasons[pid] = reason
    used_direction_ids = set()
    for d in report["directions"]:
        di = integer(d["direction_id"], "direction_id", nt*np_-1)
        if di in used_direction_ids:
            raise ValueError("Duplicate direction")
        used_direction_ids.add(di)
        if len(d["hits"]) != d["hit_count"]:
            raise ValueError("Recorded hits do not match hit_count")
        seen = set()
        causes = set()
        for h in d["hits"]:
            ci = h["compact_index"]
            if ci not in patch_by_compact or patch_by_compact[ci]["patch_id"] != h["patch_id"]:
                raise ValueError("Hit/patch witness mismatch")
            pid = h["patch_id"]
            root = integer(h["root_index"], "root_index", 1)
            if (ci, root) in seen:
                all_issues.append({"direction_id": di, "issue": "duplicate compact_index/root_index"})
            seen.add((ci, root))
            if len(h["uvtr_bits"]) != 4:
                raise ValueError("Invalid uvtr")
            u, v, t, residual = map(f32, h["uvtr_bits"])
            js = exact_by_id[pid].get("jacobians", [])
            j_sign = "unavailable"
            if js and math.isfinite(u) and math.isfinite(v):
                uu, vv = rational(h["uvtr_bits"][0]), rational(h["uvtr_bits"][1])
                weights = ((1-uu)*(1-vv), uu*(1-vv), (1-uu)*vv, uu*vv)
                j_sign = str(sign(sum(w*Fraction(q) for w,q in zip(weights,js))))
            status = patch_by_compact[ci]["stored_status"]
            if status == 3:
                causes.add(reasons[pid])
            if h["flags"] & 1:
                causes.add("hit_on_closed_boundary")
            if h["flags"] & 2:
                causes.add("singular_hit")
            hit_rows.append({"direction_id": di, "patch_id": pid, "root_index": root,
                             "stored_status": status, "hit_flags": h["flags"], "u": u, "v": v,
                             "t": t, "stored_residual": residual,
                             "exact_jacobian_sign_at_stored_uv": j_sign, "patch_reason": reasons[pid]})
        # Captured first-problem IDs may explain failures with NO accepted hit.
        for key in ("query_first_problem", "optical_first_problem", "focal_first_problem"):
            pid = d[key]
            if pid != 0xFFFFFFFF and pid in reasons:
                causes.add(reasons[pid])
        direction_rows.append({"direction_id": di, "theta_deg": math.degrees(d["theta_rad"]),
                               "phi_deg": math.degrees(d["phi_rad"]), "first_unavailable_stage": first_stage(d),
                               "query_flags": d["query_flags"], "optical_flags": d["optical_flags"],
                               "focal_flags": d["focal_flags"], "diffraction_flags": d["diffraction_flags"],
                               "hit_count": d["hit_count"], "rejected_hits": d["rejected_hits"],
                               "recorded_patch_causes": ";".join(sorted(causes))})
    if len(patch_rows) != report["summary"]["selected_patches"] or len(direction_rows) != report["summary"]["selected_directions"]:
        raise ValueError("Report summary disagrees with witness arrays")
    summary = {"format": "rainbow_patch_failure_exact_audit_v1", "origin": report["origin"],
               "captured_summary": report["summary"], "classification_counts": dict(Counter(r["classification"] for r in patch_rows)),
               "pending_patch_reason_counts": dict(Counter(r["reason"] for r in patch_rows if r["stored_status"] == 3)),
               "first_unavailable_stage_counts": dict(Counter(r["first_unavailable_stage"] for r in direction_rows)),
               "integrity_issues": all_issues, "exact_witnesses": exact_by_id,
               "absent_problem_patch_ids": report["absent_problem_patch_ids"],
               "limitations": ["Exact only for stored binary32 bilinear corners, not true optical geometry.",
                               "Host interval replays are not GPU instruction traces.",
                               "Uniform J signs plus a hemisphere certificate do not validate field/phase/area accuracy.",
                               "No counters or optical values have been modified. No repair or NaN filling is performed.",
                               "For rejected query candidates without accepted hits, only first_problem IDs were available."]}
    return summary, patch_rows, direction_rows, hit_rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        report = load(args.report)
        result, patches, directions, hits = analyze(report)
        result["source_sha256"] = hashlib.sha256(args.report.read_bytes()).hexdigest()
        args.out.mkdir(parents=True, exist_ok=False)
        with (args.out/"summary.json").open("w", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=False)
        csv_write(args.out/"patches.csv", patches,
                  ["patch_id","family","cell_x","cell_y","stored_status","classification","reason","jacobian_signs",
                   "host_orientation32","interval_orientation64","center_exactly_on_axis","interval_audit_issues"])
        csv_write(args.out/"directions.csv", directions,
                  ["direction_id","theta_deg","phi_deg","first_unavailable_stage","query_flags","optical_flags",
                   "focal_flags","diffraction_flags","hit_count","rejected_hits","recorded_patch_causes"])
        csv_write(args.out/"hits.csv", hits,
                  ["direction_id","patch_id","root_index","stored_status","hit_flags","u","v","t","stored_residual",
                   "exact_jacobian_sign_at_stored_uv","patch_reason"])
        print("Origin:", report["origin"])
        print("Selected directions:", len(directions), "; patches:", len(patches))
        print("Pending-patch reasons:", result["pending_patch_reason_counts"])
        print("First unavailable stage:", result["first_unavailable_stage_counts"])
        print("Integrity issues:", len(result["integrity_issues"]))
        print("No physical-validation PASS/FAIL or automatic repair was performed.")
        return 2 if result["integrity_issues"] else 0
    except (ValueError, KeyError, TypeError, OSError, OverflowError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
