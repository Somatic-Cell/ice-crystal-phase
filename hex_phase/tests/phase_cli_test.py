from __future__ import annotations
import hashlib
import json
from pathlib import Path
import subprocess
from process_output import diagnostic_bytes, run_checked
import sys
import tempfile
import numpy as np

exe = Path(sys.argv[1]).resolve()
repo = Path(sys.argv[2]).resolve()
sys.path.insert(0, str(repo / "python"))
from rainbow_dataset.record import PhaseRecord


def run(args: list[str], code: int = 0) -> subprocess.CompletedProcess[bytes]:
    return run_checked([str(exe), *map(str, args)], code)


def arrays(record: Path):
    with PhaseRecord(record, validate=True) as r:
        a = np.array(r.phi_cdf)
        t = np.array(r.theta_cdf)
        u = np.array(r.u_edges)
        p = np.diff(a)[:, None] * np.diff(t, axis=1)
        return a, t, u, p, r.metadata


with tempfile.TemporaryDirectory(prefix="ice_phase_検証_") as td:
    base = Path(td)
    plan = base / "orientations.csv"
    plan.write_bytes((repo / "hex_phase/configs/two_orientations.csv").read_bytes())
    common = ["--orientations", plan, "--radius-mm", "0.1", "--length-mm", "0.2", "--ki", "1", "-0.4", "0.3",
              "--ior", "1.31", "--wavelength-nm", "550", "--theta", "72", "--phi", "144", "--cdf-theta", "36", "--cdf-phi", "72"]
    run(["--help"])
    out = base / "reference"
    result = run([*common, "--out", out, "--batch-size", "512"])
    print("stdout bytes:", diagnostic_bytes(result.stdout))
    expected = {"metadata.json", "phi_cdf.npy", "theta_given_phi_cdf.npy", "u_edges.npy"}
    assert {p.name for p in out.iterdir()} == expected
    a, t, u, mass, m = arrays(out)
    assert a.shape == (73,) and t.shape == (72, 37) and u.shape == (37,)
    assert all(x.dtype.str == "<f8" and x.flags.c_contiguous for x in (a, t, u))
    assert m["orientation"]["csv_sha256"] == hashlib.sha256(plan.read_bytes()).hexdigest()
    assert m["histogram"]["relative_mass_audit_error"] < 1e-10
    assert m["sampling"]["planned_rays"] == 3072
    assert abs(mass.sum()-1) < 1e-12
    expected_g = np.sum(mass * (1-u[:-1]-u[1:])[None, :])
    assert abs(expected_g-m["hg"]["g"]) < 1e-12
    area, escaped, unresolved = (m["geometry"]["projected_area_mm2"], m["transport_audit"]["escaped_mm2"],
                               m["transport_audit"]["unresolved_mm2"])
    assert abs((escaped+unresolved)/area-1) < 1e-12
    assert m["material"]["automatic_ice_dispersion"] is False
    assert m["diffraction_policy"]["computed"] is False
    assert m["storage_processing"]["gaussian_sigma_degrees"] == 0
    assert m["quality"]["angular_convergence_certified"] is False
    # Chunking must not alter the sample IDs or orientation probabilities.
    other = base / "other_batch"
    run([*common, "--out", other, "--batch-size", "31"])
    _, _, _, other_mass, _ = arrays(other)
    error = np.abs(other_mass-mass).sum()
    assert error < 1e-10, error
    # Plateau-safe inversion, uniform in u rather than theta; independent check of sampling convention.
    rng = np.random.default_rng(1907)
    count = 100000
    up, ut = rng.random(count), rng.random(count)
    j = np.searchsorted(a, up, side="right")-1
    us = np.empty(count)
    for col in np.unique(j):
        selected = j == col
        i = np.searchsorted(t[col], ut[selected], side="right")-1
        du_cdf = t[col, i+1]-t[col, i]
        assert np.all(du_cdf > 0)
        fraction = (ut[selected]-t[col, i])/du_cdf
        us[selected] = u[i] + fraction*(u[i+1]-u[i])
    g_sample = np.mean(1-2*us)
    assert abs(g_sample-expected_g) < 6/np.sqrt(count)
    # A successfully published record is not overwritten.
    digest = {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in out.iterdir()}
    run([*common, "--out", out], code=1)
    assert digest == {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in out.iterdir()}
    for name, options in [
        ("limited", ["--max-internal-hits", "1"]),
        ("large_tail", ["--tail-tolerance", "0.001"]),
        ("coarsening", ["--max-coarsening-tv", "0"]),
        ("capacity", ["--max-output-records", "1"]),
    ]:
        failed = base / name
        run([*common, "--out", failed, *options], code=2)
        assert not failed.exists()
        part = base / (name+".part")
        assert json.loads((part / "failure.json").read_text())["complete"] is False
        assert not (part / "metadata.json").exists()
        try:
            PhaseRecord(part)
        except ValueError:
            pass
        else:
            raise AssertionError("accepted incomplete .part")
    # Matched indices: all geometric light remains in the forward direction;
    # finite-cell CDF distributes the pole mass uniformly over azimuth.
    transparent = base / "transparent"
    matched = common.copy()
    matched[matched.index("--ior")+1] = "1"
    run([*matched, "--out", transparent, "--tail-tolerance", "0"])
    pa, pt, pu, pm, mm = arrays(transparent)
    assert np.abs(pm[:, 0]-1/72).max() < 1e-13
    assert np.abs(pm[:, 1:]).max() == 0
    assert abs(mm["hg"]["g"]-(1-pu[1])) < 1e-12
    # Upstream reader bytes are not edited by this test.
    print(f"PASS CLI + unchanged PhaseRecord + NumPy: chunk_L1={error:.4g}, sampled_g={g_sample:.9g}, exact_g={expected_g:.9g}")
