#!/usr/bin/env python3
"""Generate explicit weighted SO(3) nodes. Built-ins are test/model examples,
not fits to an observed ice-crystal orientation distribution. Keep this CSV:
the record stores its SHA256, not every orientation. Standard library only.
"""
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path
import random


def multiply(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, ...]:
    w, x, y, z = a
    v, i, j, k = b
    return (w*v-x*i-y*j-z*k, w*i+x*v+y*k-z*j,
            w*j-x*k+y*v+z*i, w*k+x*j-y*i+z*v)


def axis_angle(axis: tuple[float, ...], angle: float) -> tuple[float, ...]:
    s = math.sin(angle/2)
    return (math.cos(angle/2), *(s*x for x in axis))


def quaternion(theta: float, azimuth: float, spin: float) -> tuple[float, ...]:
    # Ry(azimuth) Rz(theta) Ry(spin): the body +y axis has dot(+y)=cos(theta).
    q = multiply(multiply(axis_angle((0., 1., 0.), azimuth),
                          axis_angle((0., 0., 1.), theta)),
                 axis_angle((0., 1., 0.), spin))
    n = math.sqrt(sum(x*x for x in q))
    return tuple(x/n for x in q)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--mode", choices=("fixed", "isotropic", "axis-cone", "axis-band"), required=True)
    p.add_argument("--count", type=int, required=True)
    p.add_argument("--samples-per-orientation", type=int, required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--half-angle-deg", type=float,
                   help="cone about +y, or band about the horizontal plane; uniform solid angle")
    a = p.parse_args()
    if a.count <= 0 or a.samples_per_orientation <= 0 or a.seed < 0:
        p.error("count/samples must be positive and seed nonnegative")
    if a.samples_per_orientation >= 2**64 or a.count*a.samples_per_orientation >= 2**64:
        p.error("planned sample count exceeds uint64")
    if a.mode == "fixed" and a.count != 1:
        p.error("fixed requires count=1")
    if a.mode in ("axis-cone", "axis-band"):
        if a.half_angle_deg is None or not math.isfinite(a.half_angle_deg) or not 0 <= a.half_angle_deg <= 90:
            p.error("cone/band requires half-angle-deg in [0,90]")
    elif a.half_angle_deg is not None:
        p.error("half-angle-deg applies only to cone/band")
    pose_rng = random.Random(a.seed)
    # Separate random stream for entry-point seeds; orientation nodes do not
    # change when the ray count or ray seed implementation is modified.
    ray_rng = random.Random(a.seed ^ 0xD1B54A32D192ED03)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("x", encoding="utf-8", newline="") as f:
        f.write(f"# test/model orientation plan; mode={a.mode}; count={a.count}; seed={a.seed}; half_angle_deg={a.half_angle_deg}\n")
        f.write("# Uniform spin and azimuth; cone/band are NOT meteorological Gaussian fits.\n")
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(("qw", "qx", "qy", "qz", "weight", "samples", "seed"))
        for _ in range(a.count):
            if a.mode == "fixed":
                q = (1., 0., 0., 0.)
            else:
                u, v, w = pose_rng.random(), pose_rng.random(), pose_rng.random()
                if a.mode == "isotropic":
                    mu = 2*u-1
                elif a.mode == "axis-cone":
                    # 1-cos(alpha), stable even for a very small cone.
                    d = 2*math.sin(math.radians(a.half_angle_deg)/2)**2
                    mu = 1-u*d
                else:
                    mu = (2*u-1)*math.sin(math.radians(a.half_angle_deg))
                q = quaternion(math.acos(mu), 2*math.pi*v, 2*math.pi*w)
            writer.writerow((*[format(x, ".17g") for x in q], format(1/a.count, ".17g"),
                             a.samples_per_orientation, ray_rng.getrandbits(64)))


if __name__ == "__main__":
    main()
