"""Interpolate the statistics of the compressible impulse runs (JacobianFreeDirect) onto a points grid.

The driver writes <session>_<k>_stats.fld every rescaling period: time averages from t = 0 of the linearised
u, v, p and of the raw products uu, uv, vv, pp (not centred, unlike the incompressible ReynoldsStresses
filter). For an impulse run, average x T = int_0^T q q^T dt, so N(x) = log(A/A0) follows from the last file.

For every case this fetches mesh.xml, session.xml and two stats files (the last one and the one DK periods
earlier, as a convergence check) from the HPC with rsync, interpolates them locally onto POINTS and writes

    <local case dir>/data/pointsimp_n600.dat   x, y, u, v, p, uu, uv, vv with uu = <uu> - <u>^2 (centred), so
                                               that pp.DeltaN_computation reads it like the incompressible files
    <local case dir>/data/pointsimp_conv.json  final time, and max |N_final - N_earlier| over x

Local case dirs mirror the HPC ones under src/ (comGapRe1000Ma<M>/directLinearSolver/blowingSuction/<case>,
flatPlateRe1000Ma<M>/directLinearSolver/blowingSuction).

Usage: python scripts/compImpulseToPoints.py [--force] [Ma/case ...]   (default: every case of pp_steady_cases.txt
       on the HPC, plus the flat plates)
"""

import argparse
import contextlib
import io
import json
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from pp.DeltaN_computation import computeNx

SRC = os.path.expanduser("~/Desktop/PhD/src")
REMOTE = "hpc:Desktop/PhD/runs"
# the 2026-11-01 reruns of Ma 0.05-0.4 ran on typhoon, in ~/compImpulseRuns with the same relative paths
TYPHOON_MACH = ("0.05", "0.2", "0.4")
# cases that ran on the HPC even though their Mach number is in TYPHOON_MACH (legacy-mesh cases redone on R1,
# 2026-11-03): "Ma case" lines of runs/pp_legacyR1_cases.txt, kept locally
# and the new R1 cases of 2026-11-05 (runs/pp_newR1_cases_20261105.txt)
HPC_CASES = [os.path.expanduser("~/compImpulseStudy/legacyR1/cases23.txt"),
             os.path.expanduser("~/compImpulseStudy/legacyR1/cases_new20261105.txt")]


def host_root(rel):
    """(ssh host, runs root) holding the outputs of a case."""
    m = re.match(r"comGapRe1000Ma([0-9.]+)/", rel)
    on_hpc = set()
    for lst in HPC_CASES:
        if os.path.isfile(lst):
            on_hpc |= {f"comGapRe1000Ma{l.split()[0]}/directLinearSolver/blowingSuction/{l.split()[1]}"
                       for l in open(lst) if l.strip()}
    if m and m.group(1) in TYPHOON_MACH and rel.rstrip("/") not in on_hpc:
        return "typhoon", "compImpulseRuns"
    return "hpc", "Desktop/PhD/runs"
POINTS = os.path.expanduser("~/compImpulseStudy/points_n600.pts")
N = 600
DK = 20  # convergence check against the stats file DK rescaling periods earlier (28 time units each)
MACH = ["0.05", "0.2", "0.4", "0.6", "0.8"]


def gen_points():
    """x from -70 (upstream of the gap, for the A0 minimum) to 990 every 5, y in [0, 150] clustered at the wall."""
    if os.path.isfile(POINTS):
        return
    x = np.arange(-70, 991, 5.0)
    y = 150 * np.linspace(0, 1, N) ** 2
    X, Y = np.meshgrid(x, y, indexing="ij")
    os.makedirs(os.path.dirname(POINTS), exist_ok=True)
    np.savetxt(POINTS, np.column_stack((X.ravel(), Y.ravel())), fmt="%.6f",
               header='<?xml version="1.0" encoding="utf-8" ?>\n<NEKTAR>\n<POINTS DIM="2" FIELDS="">',
               footer="</POINTS>\n</NEKTAR>\n", comments="")


def cases():
    out = []
    lst = subprocess.run(["ssh", "-o", "BatchMode=yes", "hpc", "cat ~/Desktop/PhD/runs/pp_steady_cases.txt"],
                         capture_output=True, text=True, check=True).stdout.split("\n")
    for line in lst:
        if line.strip():
            M, c, _ = line.split()
            out.append(f"comGapRe1000Ma{M}/directLinearSolver/blowingSuction/{c}")
    out += [f"flatPlateRe1000Ma{M}/directLinearSolver/blowingSuction" for M in MACH]
    return out


def read_dat(path):
    d = np.loadtxt(path, skiprows=3)
    names = re.search(r"Variables\s*=\s*(.*)", open(path).readline()).group(1)
    return d, [v.strip() for v in names.split(",")]


def to_centred(src, dst):
    d, names = read_dat(src)
    col = {n: d[:, i] for i, n in enumerate(names)}
    u, v = col["u"], col["v"]
    out = np.column_stack((col["x"], col["y"], u, v, col["p"], col["uu"] - u * u, col["uv"] - u * v, col["vv"] - v * v))
    with open(dst, "w") as fh:
        fh.write(f"Variables = x, y, u, v, p, uu, uv, vv\n\nZone, I={len(out)}, F=POINT\n")
        np.savetxt(fh, out, fmt="%.10e")


def nx(path):
    with contextlib.redirect_stdout(io.StringIO()):
        return computeNx(path, False)


def process(rel, force=False):
    local = os.path.join(SRC, rel, "data")
    final_dat = os.path.join(local, f"pointsimp_n{N}.dat")
    if os.path.isfile(final_dat) and not force:
        return f"{rel}: exists"
    host, root = host_root(rel)
    ls = subprocess.run(["ssh", "-o", "BatchMode=yes", host, f"cd ~/{root}/{rel} && ls -d *_stats.fld 2>/dev/null"],
                        capture_output=True, text=True).stdout.split()
    ks = sorted(int(m.group(1)) for f in ls if (m := re.fullmatch(r".*_(\d+)_stats\.fld", f)))
    if not ks:
        return f"{rel}: no stats files"
    k1 = ks[-1]
    k0 = max(k for k in ks if k <= k1 - DK) if any(k <= k1 - DK for k in ks) else ks[0]
    name = re.sub(r".*/", "", rel) if "comGap" in rel else "d0_w0"  # extract_depth_width needs a d_w name
    os.makedirs(local, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        files = ["mesh.xml", "session.xml"] + [f for f in ls if re.fullmatch(rf".*_({k0}|{k1})_stats\.fld", f)]
        subprocess.run(["rsync", "-aL"] + [f"{host}:{root}/{rel}/{f}" for f in files] + [tmp + "/"], check=True)
        res = {}
        for k in (k0, k1):
            fld = [f for f in files if f.endswith(f"_{k}_stats.fld")][0]
            raw = os.path.join(tmp, f"raw_{k}.dat")
            r = subprocess.run(["FieldConvert", "-f", "-m",
                                f"interppoints:fromxml={tmp}/mesh.xml:fromfld={tmp}/{fld}:topts={POINTS}", raw],
                               cwd=tmp, capture_output=True, text=True)
            if r.returncode or not os.path.isfile(raw):
                return f"{rel}: FieldConvert failed ({r.stderr[-200:]})"
            out = os.path.join(tmp, f"{name}_{k}.dat")
            to_centred(raw, out)
            res[k] = nx(out)
            info = open(os.path.join(tmp, fld, "Info.xml")).read() if os.path.isdir(os.path.join(tmp, fld)) else ""
            if k == k1:
                os.replace(out, final_dat) if os.stat(tmp).st_dev == os.stat(local).st_dev else \
                    subprocess.run(["cp", out, final_dat], check=True)
                t_final = float(m.group(1)) if (m := re.search(r"<FinalTime>([^<]*)", info)) else float("nan")
    (x0, N0), (x1, N1) = res[k0], res[k1]
    xc = x1[(x1 >= max(x0[0], x1[0])) & (x1 <= min(x0[-1], x1[-1]))]
    dN = np.abs(np.interp(xc, x1, N1) - np.interp(xc, x0, N0))
    conv = {"k_final": k1, "k_check": k0, "t_final": t_final,
            "max_dN": float(dN.max()), "max_dN_x_le_600": float(dN[xc <= 600].max())}
    json.dump(conv, open(os.path.join(local, "pointsimp_conv.json"), "w"), indent=1)
    return f"{rel}: k {k0}->{k1}, t {t_final:.0f}, max|dN| {conv['max_dN']:.1e} (x<=600: {conv['max_dN_x_le_600']:.1e})"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("cases", nargs="*")
    a = ap.parse_args()
    gen_points()
    todo = a.cases or cases()
    with ThreadPoolExecutor(6) as ex:
        for msg in ex.map(lambda r: process(r, a.force), todo):
            print(msg, file=sys.__stdout__, flush=True)  # redirect_stdout in the worker threads is process-global
