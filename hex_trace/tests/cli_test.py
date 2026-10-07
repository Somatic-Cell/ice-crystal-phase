"""CLI contract tests; standard library only. No GPU or optical GT claimed."""
import csv
import json
import pathlib
import subprocess
import sys
import tempfile

exe = pathlib.Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory(prefix="ice_hex_cli_") as tmp:
    root = pathlib.Path(tmp)
    common = [str(exe), "--radius-mm", "0.1", "--length-mm", "0.2",
              "--ki", "0", "-1", "0", "--ior", "1.31",
              "--wavelength-nm", "550", "--samples", "17", "--batch-size", "5"]
    def invoke(extra, expected):
        p = subprocess.run(common + extra, capture_output=True, text=True)
        assert p.returncode == expected, (p.returncode, p.stdout, p.stderr)
        return p
    directory = root / "accepted"
    p = invoke(["--out", str(directory)], 0)
    m = json.loads((directory / "trace_metadata.json").read_text())
    assert m["complete"] and m["quality_accepted"] and not m["phase_cdf_record"]
    assert not m["renormalized_to_escaped_power"]
    assert m["processed_samples"] == 17
    rows = list(csv.DictReader((directory / "outgoing.csv").open()))
    assert len(rows) == m["outgoing_count"]
    recovered = sum(float(r["weight_mm2"]) for r in rows)
    expected = m["projected_area_mm2"] * m["escaped_power_sum"] / 17
    assert abs(recovered - expected) < 1e-14
    for r in rows:
        power = 0.5 * sum(float(r[k]) ** 2 for k in (
            "j00_re", "j00_im", "j10_re", "j10_im", "j01_re", "j01_im", "j11_re", "j11_im"))
        assert abs(power - float(r["power_fraction"])) < 2e-14
    invoke(["--out", str(directory)], 1)  # no overwrite
    bad = root / "rejected"
    invoke(["--out", str(bad), "--max-internal-hits", "1"], 2)
    assert not bad.exists()
    bad_meta = json.loads((root / "rejected.part" / "trace_metadata.json").read_text())
    assert not bad_meta["complete"] and not bad_meta["quality_accepted"]
    assert bad_meta["status_counts"]["interaction_limit"] > 0
    invoke(["--out", str(root / "invalid"), "--tail-tolerance", "1"], 1)
    invoke(["--out", str(root / "cuda"), "--backend", "cuda"], 1)  # CPU exe never falls back
    invoke(["--out", str(root / "duplicate"), "--ior", "1.31"], 1)
    assert not (root / "invalid.part").exists()
    print("PASS CLI: accepted output, CSV Jones/area weights, refused overwrite, rejected tail, invalid/duplicate input, no CUDA fallback")
    print(p.stdout.strip())
