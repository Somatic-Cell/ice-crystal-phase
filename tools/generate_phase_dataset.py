#!/usr/bin/env python3
"""Sequential, resumable condition-grid generation; no point sampling or learning.

Standard library only. Every condition invokes the real rainbow_generate once.
This orchestrator does not batch OptiX queries or share GPU contexts across jobs.
It never changes resolution/tolerances to get a failed condition to pass.
"""
from __future__ import annotations

import argparse
import ast
from contextlib import contextmanager
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import time
from typing import Any, Iterator

ROOT = Path(__file__).resolve().parents[1]
VERSION = "rainbow.condition_batch.v1"
GENERATOR_VERSION = "numpy_cdf_water_geometry_v1"
MODULES = ("raindrop_trace.optixir", "patch_build.fatbin", "patch_query.optixir",
           "patch_optics.fatbin", "phase_cdf.fatbin")
ARRAYS = ("phi_cdf.npy", "theta_given_phi_cdf.npy", "u_edges.npy")


def fail_constant(value: str) -> None:
    raise ValueError(f"Nonfinite JSON constant: {value}")


def load_json(path: Path, limit: int = 8 << 20) -> dict[str, Any]:
    if path.stat().st_size > limit:
        raise ValueError(f"JSON too large: {path}")
    with path.open("r", encoding="utf-8-sig") as stream:
        value = json.load(stream, parse_constant=fail_constant)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def exact(value: Any) -> Decimal:
    if isinstance(value, bool):
        raise ValueError("Boolean is not a numerical parameter")
    try:
        d = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError(f"Invalid decimal: {value}") from exc
    if not d.is_finite():
        raise ValueError("Nonfinite parameter")
    return d


def scalar(value: Any, low: float, high: float) -> float:
    d = exact(value)
    if not Decimal(str(low)) <= d <= Decimal(str(high)):
        raise ValueError(f"Parameter {value} is outside [{low}, {high}]")
    return float(d)


def integer(value: Any, low: int, high: int) -> int:
    d = exact(value)
    if d != d.to_integral_value() or not low <= d <= high:
        raise ValueError(f"Expected integer in [{low}, {high}]")
    return int(d)


def flag(value: Any) -> bool:
    if type(value) is not bool:
        raise ValueError("Expected JSON true or false")
    return value


def decimal_string(d: Decimal) -> str:
    return format(d, "f")


def wavelengths(spec: dict[str, Any]) -> list[str]:
    if set(spec) != {"start_nm", "stop_nm", "step_nm"}:
        raise ValueError("wavelengths needs start_nm, stop_nm, step_nm")
    start, stop, step = (exact(spec[k]) for k in ("start_nm", "stop_nm", "step_nm"))
    if not 380 <= start <= stop <= 830 or step <= 0:
        raise ValueError("Require 380 <= start <= stop <= 830 and positive step")
    intervals = (stop-start)/step
    if intervals != intervals.to_integral_value():
        raise ValueError("The wavelength range must contain an integer number of steps; no silent endpoint truncation")
    count = int(intervals)+1
    if count > 100000:
        raise ValueError("Too many wavelengths")
    values = [decimal_string(start+step*i) for i in range(count)]
    packed = [struct.pack("<f", float(v)) for v in values]
    if len(set(packed)) != len(packed):
        raise ValueError("Distinct wavelengths collapse to the same FP32 solver value")
    return values


def inclinations(spec: dict[str, Any]) -> list[str]:
    if set(spec) != {"min_degrees", "max_degrees", "count", "placement"}:
        raise ValueError("inclinations needs min_degrees, max_degrees, count, placement")
    if spec["placement"] != "linear_angle_cell_centers":
        raise ValueError("Only explicit linear_angle_cell_centers is supported")
    low, high = exact(spec["min_degrees"]), exact(spec["max_degrees"])
    count = integer(spec["count"], 1, 100000)
    if not -90 <= low < high <= 90:
        raise ValueError("Canonical inclination domain is [-90,90], NOT [0,180]")
    return [decimal_string(low+(high-low)*(Decimal(i)+Decimal("0.5"))/count)
            for i in range(count)]


def normalized_config(raw: dict[str, Any]) -> dict[str, Any]:
    allowed = {"output_root", "generator", "radius_mm", "force_sphere", "temperature_c",
               "pressure_pa", "wavelengths", "inclinations", "solver", "reserve_gib"}
    if set(raw) != allowed:
        raise ValueError(f"Unexpected/missing config keys: {set(raw)^allowed}")
    s = raw["solver"]
    sk = {"grid", "query_theta", "query_phi", "cdf_theta", "cdf_phi", "stage", "device",
          "focal_offsets", "allow_underresolved", "maximum_coarsening_tv",
          "cdf_max_l1", "cdf_max_lost_mass"}
    if not isinstance(s, dict) or set(s) != sk:
        raise ValueError("Invalid solver keys")
    out = {k: integer(s[k], 2 if k == "grid" else 1, 32767)
           for k in ("grid", "query_theta", "query_phi", "cdf_theta", "cdf_phi")}
    if out["grid"]**2*4 > 0xFFFFFFFF:
        raise ValueError("Source index range exceeds the existing solver contract")
    for axis in ("theta", "phi"):
        if out["query_"+axis] % out["cdf_"+axis] or out["cdf_"+axis] > out["query_"+axis]:
            raise ValueError("CDF grids must divide query grids")
    if s["stage"] not in ("diffraction", "focal", "path", "incoherent"):
        raise ValueError("Unknown stage")
    out["stage"] = s["stage"]
    out["device"] = integer(s["device"], 0, 128)
    out["allow_underresolved"] = flag(s["allow_underresolved"])
    out["maximum_coarsening_tv"] = scalar(s["maximum_coarsening_tv"], 0, 1)
    out["cdf_max_l1"] = scalar(s["cdf_max_l1"], 0, 1)
    out["cdf_max_lost_mass"] = scalar(s["cdf_max_lost_mass"], 0, 1)
    offsets = s["focal_offsets"]
    if not isinstance(offsets, list) or len(offsets) != 4:
        raise ValueError("Four focal offsets required")
    out["focal_offsets"] = [integer(v, 0, 3) for v in offsets]
    radius = scalar(raw["radius_mm"], 1e-9, 3)
    if s["stage"] in ("focal", "diffraction") and not .1 <= radius <= 1:
        raise ValueError("Table-II diffraction path currently accepts radii 0.1..1.0 mm only")
    return {"radius_mm": radius, "force_sphere": flag(raw["force_sphere"]),
            "temperature_c": scalar(raw["temperature_c"], 0, 60),
            "pressure_pa": scalar(raw["pressure_pa"], 50000, 100000000),
            "wavelength_nm": wavelengths(raw["wavelengths"]),
            "inclination_degrees": inclinations(raw["inclinations"]),
            "inclination_placement": "linear_angle_cell_centers",
            "material": "water", "solver": out}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def object_hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode("utf-8")).hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix+".tmp")
    with tmp.open("w", encoding="utf-8", newline="\n") as f:
        json.dump(value, f, indent=2, ensure_ascii=False, allow_nan=False)
        f.write("\n")
    tmp.replace(path)


@contextmanager
def batch_lock(root: Path) -> Iterator[None]:
    lock = root / ".batch.lock"
    try:
        with lock.open("x", encoding="ascii") as f:
            import os
            f.write(f"pid={os.getpid()}\n")
    except FileExistsError as exc:
        raise RuntimeError("Batch lock exists. Do not run two writers. After verifying that no batch is running, remove the stale .batch.lock manually.") from exc
    try:
        yield
    finally:
        lock.unlink()


def npy_size(path: Path, expected_shape: tuple[int, ...]) -> int:
    """Header and length validation only; no full CDF read, no sampling."""
    with path.open("rb") as f:
        if f.read(6) != b"\x93NUMPY":
            raise ValueError(f"Not NPY: {path}")
        ver = f.read(2)
        if ver == b"\x01\x00":
            raw = f.read(2)
            if len(raw) != 2: raise ValueError("Truncated NPY")
            length = struct.unpack("<H", raw)[0]
        elif ver == b"\x02\x00":
            raw = f.read(4)
            if len(raw) != 4: raise ValueError("Truncated NPY")
            length = struct.unpack("<I", raw)[0]
        else:
            raise ValueError("Unsupported NPY version")
        if not 0 < length <= 65536: raise ValueError("Invalid NPY header size")
        text = f.read(length)
        if len(text) != length: raise ValueError("Truncated NPY header")
        h = ast.literal_eval(text.decode("latin1"))
        if (h.get("descr") != "<f8" or h.get("fortran_order") is not False
                or h.get("shape") != expected_shape):
            raise ValueError(f"Wrong NPY layout: {path}")
        expected = f.tell()+8*math.prod(expected_shape)
    if path.stat().st_size != expected:
        raise ValueError(f"Truncated or trailing NPY data: {path}")
    return expected


def close(a: Any, b: float, tol: float = 1e-12) -> bool:
    return (type(a) in (int, float) and math.isfinite(a)
            and abs(a-b) <= tol*max(1.0, abs(b)))


def fp32(x: float) -> float:
    return struct.unpack("<f", struct.pack("<f", x))[0]


def verify_record(path: Path, plan: dict[str, Any], ii: int, ll: int) -> dict[str, Any]:
    c, s = plan["conditions"], plan["conditions"]["solver"]
    m = load_json(path/"metadata.json", 1 << 20)
    def require(ok: bool, message: str) -> None:
        if not ok: raise ValueError(f"{path}: {message}")
    expected = {"complete": True, "schema": "rainbow.phase_cdf.numpy.v2",
                "generator_version": GENERATOR_VERSION,
                "coordinate_contract_id": "rainbow.phase_cdf.coordinates.v1",
                "direction_convention": "physical_propagation", "dtype": "<f8", "order": "C",
                "density_measure": "solid_angle_sr", "cdf_axis_order": ["phi_cell", "theta_edge"],
                "theta_count": s["cdf_theta"], "phi_count": s["cdf_phi"],
                "query_grid": [s["query_theta"], s["query_phi"]],
                "incident_grid": [s["grid"], s["grid"]], "stage": s["stage"],
                "force_sphere": c["force_sphere"], "focal_quarter_turn_offsets": s["focal_offsets"]}
    for key, value in expected.items(): require(m.get(key) == value, f"mismatched {key}")
    require(close(m.get("radius_mm"), fp32(c["radius_mm"])), "radius")
    require(close(m.get("wavelength_nm"), fp32(float(c["wavelength_nm"][ll]))), "wavelength")
    alpha = float(c["inclination_degrees"][ii])
    require(close(m.get("incident_inclination_degrees"), alpha), "inclination")
    material = m.get("material", {})
    require(material.get("model") == "water_iapws_r9_97", "material")
    require(close(material.get("temperature_kelvin"), c["temperature_c"]+273.15), "temperature")
    require(close(material.get("pressure_pascal"), c["pressure_pa"]), "pressure")
    mat = m.get("sampling_frame_columns", [])
    require(isinstance(mat, list) and len(mat)==3 and all(isinstance(row,list) and len(row)==3 for row in mat), "frame shape")
    require(all(type(v) in (int,float) and math.isfinite(v) for row in mat for v in row), "frame values")
    cols = [[mat[r][j] for r in range(3)] for j in range(3)]
    for j in range(3):
        for k in range(3):
            require(close(sum(cols[j][r]*cols[k][r] for r in range(3)), float(j==k), 2e-12), "frame orthogonality")
    cross = [cols[0][1]*cols[1][2]-cols[0][2]*cols[1][1],
             cols[0][2]*cols[1][0]-cols[0][0]*cols[1][2],
             cols[0][0]*cols[1][1]-cols[0][1]*cols[1][0]]
    require(sum(cross[r]*cols[2][r] for r in range(3)) > 1-2e-12, "frame handedness")
    angle = math.radians(alpha)
    ki = [math.cos(angle), -math.sin(angle), 0.0]
    require(all(abs(cols[2][r]-ki[r]) < 2e-6 for r in range(3)), "incident direction/sign")
    geometry, cs = m.get("geometry", {}), m.get("cross_sections", {})
    area = geometry.get("projected_area_mm2")
    require(type(area) in (int,float) and math.isfinite(area) and area>0, "projected area")
    require(geometry.get("projected_area_converged") is True and geometry.get("convexity_verified") is True, "area validation")
    require(cs.get("model")=="geometric_nonabsorbing" and cs.get("wave_scattering_cross_section_computed") is False, "cross-section model")
    require(cs.get("scattering_mm2")==area and cs.get("extinction_mm2")==area and cs.get("absorption_mm2")==0, "cross sections")
    q = m.get("quality", {})
    require(q.get("invalid_values")==0 and q.get("incomplete_directions")==0, "invalid optical data")
    require(q.get("allow_underresolved") is s["allow_underresolved"], "underresolution policy")
    require(s["allow_underresolved"] or q.get("underresolved_directions")==0, "underresolution not authorized")
    for key, cfg in (("cdf_l1_mass_error","cdf_max_l1"),("cdf_lost_probability_mass","cdf_max_lost_mass")):
        value = q.get(key)
        require(type(value) in (int,float) and math.isfinite(value) and 0<=value<=s[cfg], key)
    sp = m.get("storage_processing", {})
    require(sp.get("additional_filter")=="none", "unexpected global smoothing")
    tv = sp.get("coarsening_tv")
    require(type(tv) in (int,float) and math.isfinite(tv) and 0<=tv<=s["maximum_coarsening_tv"], "coarsening error")
    hg = m.get("hg", {})
    require(hg.get("target")=="saved_cdf" and type(hg.get("g")) in (float,int) and -1<hg["g"]<1, "HG target")
    nt, np_ = s["cdf_theta"], s["cdf_phi"]
    sizes = {name: npy_size(path/name, shape) for name,shape in
             zip(ARRAYS, ((np_+1,), (np_,nt+1), (nt+1,)))}
    return {"metadata_sha256": sha256(path/"metadata.json"), "npy_sizes": sizes,
            "hg_g": hg["g"], "projected_area_mm2": area,
            "underresolved_directions": q.get("underresolved_directions"), "coarsening_tv": tv}


def command(exe: Path, out: Path, c: dict[str, Any], ii: int, ll: int) -> list[str]:
    s = c["solver"]
    values = {"out": str(out), "material": "water", "radius-mm": str(c["radius_mm"]),
              "wavelength-nm": c["wavelength_nm"][ll], "inclination-deg": c["inclination_degrees"][ii],
              "temperature-c": str(c["temperature_c"]), "pressure-pa": str(c["pressure_pa"]),
              "grid": s["grid"], "query-theta": s["query_theta"], "query-phi": s["query_phi"],
              "cdf-theta": s["cdf_theta"], "cdf-phi": s["cdf_phi"], "stage": s["stage"],
              "device": s["device"], "max-coarsening-tv": s["maximum_coarsening_tv"],
              "cdf-max-l1": s["cdf_max_l1"], "cdf-max-lost-mass": s["cdf_max_lost_mass"],
              "focal-offsets": ",".join(map(str,s["focal_offsets"]))}
    args = [str(exe)]
    for k,v in values.items(): args += ["--"+k,str(v)]
    if c["force_sphere"]: args.append("--sphere")
    if s["allow_underresolved"]: args.append("--allow-underresolved")
    return args


def run(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--out-root", type=Path)
    ap.add_argument("--generator", type=Path)
    ap.add_argument("--wavelength-range", nargs=3, metavar=("START_NM","STOP_NM","STEP_NM"))
    ap.add_argument("--radius-mm")
    ap.add_argument("--inclination-count", type=int)
    ap.add_argument("--allow-underresolved", action="store_true", help="explicitly permit recorded underresolution warnings")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--limit", type=int, help="run at most this many NEW records, without changing the planned dataset")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--execute", action="store_true", help="actually start the simulations; default is a dry run")
    mode.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    if args.limit is not None and args.limit<=0: raise ValueError("limit must be positive")
    raw = load_json(args.config)
    if args.wavelength_range: raw["wavelengths"] = dict(zip(("start_nm","stop_nm","step_nm"),args.wavelength_range))
    if args.radius_mm is not None: raw["radius_mm"]=args.radius_mm
    if args.inclination_count is not None: raw["inclinations"]["count"]=args.inclination_count
    if args.allow_underresolved: raw["solver"]["allow_underresolved"]=True
    c = normalized_config(raw)
    root = (args.out_root or Path(raw["output_root"]))
    exe = (args.generator or Path(raw["generator"]))
    root = (root if root.is_absolute() else ROOT/root).resolve()
    exe = (exe if exe.is_absolute() else ROOT/exe).resolve()
    reserve = int(scalar(raw["reserve_gib"], 0, 1000000)*(1<<30))
    nt,np_ = c["solver"]["cdf_theta"],c["solver"]["cdf_phi"]
    payload = 8*((np_+1)+np_*(nt+1)+(nt+1))
    estimate_per_record = payload+(128<<10) # allowance for headers/metadata/receipt; logs extra
    ni,nl = len(c["inclination_degrees"]),len(c["wavelength_nm"])
    print(f"conditions={ni} x {nl} = {ni*nl}; CDF={nt}x{np_}; array_bytes_per_record={payload}")
    print(f"array_total_TB={payload*ni*nl/1e12:.6f}; output={root}")
    print(f"inclination_centers=[{c['inclination_degrees'][0]}, {c['inclination_degrees'][-1]}] degrees")
    print(f"wavelengths=[{c['wavelength_nm'][0]}, {c['wavelength_nm'][-1]}] nm; inclusive exact-step endpoints")
    print("One subprocess per condition; one GPU job at a time. No context reuse / query batching added.")
    if c["solver"]["allow_underresolved"]:
        print("WARNING: underresolved diffraction is explicitly permitted, not certified accurate.")
    if c["solver"]["maximum_coarsening_tv"]==1:
        print("WARNING: maximum_coarsening_tv=1 is not an accuracy acceptance threshold.")
    if not args.execute:
        print("DRY RUN: no simulations or files written. Add --execute to start.")
        print(subprocess.list2cmdline(command(exe,root/"records/i0000/w0000",c,0,0)))
        return 0
    binaries = {"generator": exe, **{name: exe.parent/"modules"/name for name in MODULES}}
    for name,path in binaries.items():
        if not path.is_file(): raise FileNotFoundError(f"Missing {name}: {path}")
    spec = {"version": VERSION, "conditions": c,
            "binary_sha256": {name: sha256(path) for name,path in binaries.items()}}
    signature = object_hash(spec)
    plan = {**spec, "signature": signature, "record_layout": "records/i{incident_index:04d}/w{wavelength_index:04d}"}
    root.mkdir(parents=True,exist_ok=True)
    with batch_lock(root):
        manifest = root/"batch_plan.json"
        if manifest.exists():
            if not args.resume: raise FileExistsError("Dataset plan exists; use --resume with identical configuration and binaries")
            if load_json(manifest)!=plan: raise ValueError("Resume rejected: configuration/binaries differ from batch_plan.json")
        else:
            if any(p.name!=".batch.lock" for p in root.iterdir()):
                raise ValueError("Refusing to adopt a nonempty directory without this batch plan")
            write_json(manifest,plan)
        pending=[]; skipped=0
        for ii in range(ni):
            for ll in range(nl):
                job=f"i{ii:04d}_w{ll:04d}";path=root/f"records/i{ii:04d}/w{ll:04d}"
                receipt=root/"receipts"/(job+".json")
                if Path(str(path)+".part").exists():
                    raise ValueError(f"Incomplete .part exists: {path}.part; inspect and move it manually before retrying")
                if path.exists():
                    if not args.resume or not receipt.exists():
                        raise ValueError(f"Existing record without authorized resume/receipt: {path}")
                    saved=load_json(receipt);actual=verify_record(path,plan,ii,ll)
                    if saved.get("signature")!=signature or saved.get("incident_index")!=ii or saved.get("wavelength_index")!=ll or saved.get("record")!=actual:
                        raise ValueError(f"Resume receipt/content mismatch: {path}")
                    skipped+=1
                else:
                    if receipt.exists(): raise ValueError(f"Receipt exists but record missing: {path}")
                    pending.append((ii,ll,job,path,receipt))
        selected=pending if args.limit is None else pending[:args.limit]
        needed=len(selected)*estimate_per_record+reserve
        free=shutil.disk_usage(root).free
        print(f"verified_existing={skipped}; pending={len(pending)}; this_run={len(selected)}; free_TB={free/1e12:.3f}")
        if needed>free: raise RuntimeError(f"Insufficient free space: estimated {needed} bytes including reserve; available {free}")
        start=time.monotonic()
        for n,(ii,ll,job,path,receipt) in enumerate(selected,1):
            free=shutil.disk_usage(root).free
            if free<estimate_per_record+reserve: raise RuntimeError("Free-space reserve reached; stopping without lowering precision")
            (root/"logs").mkdir(exist_ok=True)
            log=root/"logs"/(job+".log")
            if log.exists():
                log=root/"logs"/(job+f".{time.time_ns()}.log")
            path.parent.mkdir(parents=True,exist_ok=True)
            cmd=command(exe,path,c,ii,ll)
            print(f"[{n}/{len(selected)}] {job}: alpha={c['inclination_degrees'][ii]} lambda={c['wavelength_nm'][ll]}",flush=True)
            with log.open("xb") as stream:
                rc=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT,check=False).returncode
            if rc!=0:
                raise RuntimeError(f"Generator exit={rc}; stopped at {job}. See {log}. No failed record was skipped or deleted.")
            checked=verify_record(path,plan,ii,ll)
            write_json(receipt,{"signature":signature,"incident_index":ii,"wavelength_index":ll,"record":checked})
            write_json(root/"progress.json",{"completed":skipped+n,"total":ni*nl,"last_record":job,
                                            "elapsed_seconds_this_run":time.monotonic()-start})
            print(f"  saved; g={checked['hg_g']:.9g}; Aproj_mm2={checked['projected_area_mm2']:.9g}",flush=True)
        complete=skipped+len(selected)==ni*nl
        write_json(root/"progress.json",{"completed":skipped+len(selected),"total":ni*nl,"complete":complete,
                                        "elapsed_seconds_this_run":time.monotonic()-start})
        print("Dataset complete." if complete else "Requested subset finished. Repeat with --resume to continue the same plan.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(run())
    except KeyboardInterrupt:
        print("Interrupted. No automatic cleanup of records/.part directories.",file=sys.stderr)
        raise SystemExit(130)
    except (OSError, ValueError, RuntimeError, KeyError, TypeError, struct.error) as exc:
        print(f"generate_phase_dataset: {exc}",file=sys.stderr)
        raise SystemExit(1)
