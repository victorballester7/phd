"""
Collect the single-mode Arnoldi beta scan with b = 8d (beta_k = 2 pi k / b, k = 1..31) of
runs/quasi3DGapRe1000inc/arnoldiLinearSolver/singleMode/betaScan_b8d/<case>/k<NN>/mesh.evl
(read through the sshfs mount of the HPC home) into data/stability/betaScan_b8d.npz.

For each run the last Arnoldi iteration of mesh.evl is used. The leading eigenvalue is the one
of largest growth among the Ritz values with residual < RES_TOL; real Ritz values with |mu| = 1 to
round-off are the spurious constant mode and are skipped.
If no Ritz value is converged, the largest-growth one is kept and flagged (converged = False).

Also collects the old full-domain d4_w15 single-mode scan (tau = 1750, arbitrary lambda_z) of
singleMode/d4_w15/Lz<x>d/NumSteps250e3 for comparison.

Usage:  cd scripts/pp && .venv/bin/python scripts/betaScanB8dCollect.py
"""

import glob
import os
import re

import numpy as np

# HPC sshfs mount by default; override with BETASCAN_RUNS=<local rsync mirror of singleMode/>
RUNS = os.path.expanduser(os.environ.get(
    "BETASCAN_RUNS", "~/hosts/hpc/PhD/runs/quasi3DGapRe1000inc/arnoldiLinearSolver/singleMode"))
SCAN = os.path.join(RUNS, "betaScan_b8d")
OUT = os.path.expanduser("~/Desktop/PhD/data/stability/betaScan_b8d.npz")
CASES = {"d4_w15": (4.0, 15.0), "d3_w21": (3.0, 21.0), "d1.75_w40": (1.75, 40.0)}
K = np.arange(1, 32)
RES_TOL = 5e-4
RES_TOL_ALL = 1e-3  # looser tolerance for the sub-dominant modes stored in <case>_modes


def leading(evl):
    """(growth, frequency, residual, iteration, converged) of the leading Ritz value, or None."""
    if not os.path.isfile(evl):
        return None
    blocks = open(evl).read().split("-- Iteration")
    last = blocks[-1]
    it = re.match(r"\s*=\s*(\d+)", last)
    rows = []
    for line in last.splitlines():
        if line.startswith("EV:"):
            p = line.split()
            mag, ang, g, f, r = map(float, p[2:7])
            if abs(mag - 1) < 1e-8 and ang == 0:  # spurious |mu| = 1 (constant) mode
                continue
            rows.append((g, f, r))
    if not rows:
        return None
    conv = [r for r in rows if r[2] < RES_TOL]
    g, f, r = max(conv or rows, key=lambda v: v[0])
    return g, abs(f), r, int(it.group(1)) if it else -1, bool(conv)


def modes(evl):
    """All Ritz values (growth, frequency >= 0, residual) of the last iteration with residual < RES_TOL_ALL,
    one per complex-conjugate pair, spurious |mu| = 1 mode skipped."""
    if not os.path.isfile(evl):
        return []
    rows = []
    for line in open(evl).read().split("-- Iteration")[-1].splitlines():
        if line.startswith("EV:"):
            mag, ang, g, f, r = map(float, line.split()[2:7])
            if (abs(mag - 1) < 1e-8 and ang == 0) or r >= RES_TOL_ALL or f < 0:
                continue
            rows.append((g, f, r))
    return rows


def main():
    out = {"k": K}
    for case, (d, w) in CASES.items():
        res = np.full((len(K), 4), np.nan)
        conv = np.zeros(len(K), bool)
        for i, k in enumerate(K):
            r = leading(os.path.join(SCAN, case, f"k{k:02d}", "mesh.evl"))
            if r:
                res[i] = r[:4]
                conv[i] = r[4]
        out[f"{case}_sigma"], out[f"{case}_omega"] = res[:, 0], res[:, 1]
        out[f"{case}_res"], out[f"{case}_iter"], out[f"{case}_conv"] = res[:, 2], res[:, 3], conv
        out[f"{case}_beta"] = 2 * np.pi * K / (8 * d)
        # all converged modes: columns k, sigma, omega, residual
        allm = [(k, *m) for k in K for m in modes(os.path.join(SCAN, case, f"k{k:02d}", "mesh.evl"))]
        out[f"{case}_modes"] = np.array(allm).reshape(-1, 4)
        n = np.isfinite(res[:, 0]).sum()
        print(f"{case}: {n}/31 runs with eigenvalues, {conv.sum()} converged", end="")
        if n:
            j = np.nanargmax(res[:, 0])
            print(f"; max sigma = {res[j, 0]:.4e} at k = {K[j]} (beta = {2 * np.pi * K[j] / (8 * d):.3f}, "
                  f"lambda_z = {8 / K[j]:.3f}d, omega = {res[j, 1]:.4f})")
        else:
            print()
    # old full-domain d4_w15 scan (tau = 1750)
    old = []
    for lz in glob.glob(os.path.join(RUNS, "d4_w15", "Lz*d")):
        r = leading(os.path.join(lz, "NumSteps250e3", "mesh.evl"))
        if r and r[4]:
            old.append((float(os.path.basename(lz)[2:-1]), r[0]))
    old = np.array(sorted(old)).reshape(-1, 2)
    out["old_d4_w15_lzd"], out["old_d4_w15_sigma"] = old[:, 0], old[:, 1]
    # validation: k = 16 on the full-length (outflow 1000) mesh
    v = leading(os.path.join(SCAN, "d4_w15_L1000_k16", "mesh.evl"))
    out["val_L1000_k16"] = np.array(v[:4] if v else [np.nan] * 4)
    if v:
        print(f"validation d4_w15 k = 16, full-length mesh: sigma = {v[0]:.4e} (converged {v[4]}); "
              f"old full-domain lambda_z = 0.5d: {dict(old.tolist()).get(0.5, np.nan):.4e}")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    np.savez(OUT, **out)
    print("->", OUT)


if __name__ == "__main__":
    main()
