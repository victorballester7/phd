"""Sum the per-node impulse runs into the T -> infinity statistics of the awgn noise run.

With awgn drawn independently at every quadrature node, the stationary covariance of the noise run
is C = sigma^2 sum_{p,c} int q_pc q_pc^T dt (``scripts/makeImpulseSourceCases.py``): the raw second
moments of the per-node impulse runs add. Each run's statistics file holds the time average from
t = 0 over the same T, i.e. mean and variance, so its raw moments are

    <u u> = U^2 + uu,   <u v> = U V + uv,   <v v> = V^2 + vv,

and the sum is written as a file with zero mean and these moments in uu, uv, vv, readable by
``computeNx`` like any noise run:

    <parent>_srcsum/data/pointsavg_n600.dat

Usage: python scripts/sumImpulseSources.py   (after scripts/impulseStatsToPoints.py)
"""

import os
import re

import numpy as np

SRC = os.path.expanduser("~/Desktop/PhD/src")
BASES = [
    os.path.join(SRC, "incGapRe1000/directLinearSolver/blowingSuction"),
    os.path.join(SRC, "flatPlateRe1000inc/directLinearSolver"),
]
N = 600


def main():
    for base in BASES:
        groups = {}
        for d in sorted(os.listdir(base)):
            m = re.fullmatch(r"(.+)_src\d+[uv]", d)
            if m:
                groups.setdefault(m.group(1), []).append(d)
        for parent, runs in groups.items():
            files = [os.path.join(base, r, "data", f"pointsavg_n{N}.dat") for r in runs]
            missing = [r for r, f in zip(runs, files) if not os.path.isfile(f)]
            if missing:
                print(f"{parent}: missing {missing}, skipped")
                continue
            total = None
            for f in files:
                d = np.loadtxt(f, skiprows=3)  # x, y, u, v, p, uu, uv, vv
                raw = np.stack([d[:, 2] ** 2 + d[:, 5], d[:, 2] * d[:, 3] + d[:, 6], d[:, 3] ** 2 + d[:, 7]], axis=1)
                if total is None:
                    xy, total = d[:, :2], raw
                else:
                    assert np.array_equal(xy, d[:, :2]), f
                    total += raw
            out = np.column_stack([xy, np.zeros((len(xy), 3)), total])
            dst = os.path.join(base, f"{parent}_srcsum", "data")
            os.makedirs(dst, exist_ok=True)
            with open(files[0]) as fh:
                header = "".join(fh.readline() for _ in range(3))
            path = os.path.join(dst, f"pointsavg_n{N}.dat")
            with open(path, "w") as fh:
                fh.write(header)
                np.savetxt(fh, out, fmt="%.8e")
            print(f"{parent}: summed {len(runs)} node impulses -> {path}")


if __name__ == "__main__":
    main()
