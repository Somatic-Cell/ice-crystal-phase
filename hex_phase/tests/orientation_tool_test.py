from __future__ import annotations
import csv
import io
import math
from pathlib import Path
import subprocess
import sys
import tempfile

script = Path(sys.argv[1]).resolve()

def load(p):
    return list(csv.DictReader(io.StringIO("\n".join(s for s in p.read_text().splitlines() if not s.startswith("#")))))

def invoke(args, code=0):
    p = subprocess.run([sys.executable, str(script), *map(str, args)], capture_output=True, text=True)
    assert p.returncode == code, p.stderr

with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    for mode in ("fixed", "isotropic", "axis-cone", "axis-band"):
        count = 1 if mode == "fixed" else 20000
        args = ["--mode", mode, "--count", count, "--samples-per-orientation", "17", "--seed", "23"]
        if mode in ("axis-cone", "axis-band"):
            args += ["--half-angle-deg", "3"]
        a, b = td/(mode+".csv"), td/(mode+"_copy.csv")
        invoke([*args, "--out", a]); invoke([*args, "--out", b]); assert a.read_bytes() == b.read_bytes()
        rows = load(a)
        assert len(rows) == count
        assert abs(math.fsum(float(r["weight"]) for r in rows)-1) < 1e-14
        mus, xx = [], []
        for r in rows:
            w, x, y, z = (float(r[k]) for k in ("qw", "qx", "qy", "qz"))
            assert abs(w*w+x*x+y*y+z*z-1) < 1e-14
            mus.append(1-2*(x*x+z*z))  # world y component of rotated body +y
            xx.append(1-2*(y*y+z*z))
            assert int(r["samples"]) == 17 and 0 <= int(r["seed"]) < 2**64
        if mode == "isotropic":
            assert abs(math.fsum(mus)/count) < .02
            assert abs(math.fsum(x*x for x in mus)/count-1/3) < .02
            assert abs(math.fsum(xx)/count) < .02  # spin is not fixed
        elif mode == "axis-cone":
            lo = math.cos(math.radians(3))
            assert min(mus) >= lo-1e-14 and max(mus) <= 1+1e-14
            assert abs(math.fsum(mus)/count-(1+lo)/2) < (1-lo)*.02
        elif mode == "axis-band":
            hi = math.sin(math.radians(3))
            assert max(map(abs, mus)) <= hi+1e-14
            assert abs(math.fsum(mus)/count) < hi*.02
    print("PASS orientation plan generator: four modes, reproducibility, unit quaternions, solid-angle and spin moments")
