"""
Local (parallel flow) spatial stability of three gaps shared by every Mach number, d3_w13, d1.5_w25 and
d0.75_w90 (Re_delta*_le = 1000), at Ma = 0.05 ... 0.8 and incompressible: the frozen base-flow profile at every
station x is solved for a band of real reduced frequencies F, as in scripts/predictNfactorLocal.py, but with the
compressible operator (pp.compressibleStability) and the thermodynamics of the Nektar++ sessions (gamma = 1.4,
Pr = 0.72, Sutherland S = 110 K / Tref, Tref = 225 K, isothermal wall T' = 0). Incompressible = Ma 0 with constant
viscosity (Orr-Sommerfeld exactly). Mode selection as pp.localNfactor.select_mode: least damped mode with phase speed
in (0.15, 0.75), |alpha_i| < 0.9 alpha_r.

Base flows (steady DNS checkpoints):
    inc            runs/incGapRe1000/baseflow/dns/<case>  (d3_w13 -> d3_w13_randIC), converged to < 2e-8
    Ma0.05-0.4,    runs/comGapRe1000Ma<M>/baseflow/dns/<case>, the checkpoint of runs/pp_steady_cases.txt
    Ma0.8 d1.5_w25, d0.75_w90
    Ma0.6 (all),   LEGACY runs (NewtonRelativeIteTol 1e-10, domain W+500 for d3_w13/d1.5_w25), last checkpoint,
    Ma0.8 d3_w13   not converged: velocity swing over the last 150 tu <= 5e-3 in the gap; downstream Ma0.6
                   d1.5_w25 exceeds 1e-2 from x ~ 300 and Ma0.6 d0.75_w90 carries a TS packet of amplitude up to
                   0.3 at x ~ 400-600, so their stations stop at XMAX_OVERRIDE (swing <= 1e-2 there).

Stages
    --extract   rsync mesh + checkpoint, FieldConvert interppoints on the station lines -> data/lst/sections/*.npz
    --solve     alpha(x, F) per (flow, case) -> data/lst/gap_<flow>_<case>.npz (x, F, alpha)
usage: python scripts/lstGapMach.py --extract --solve [--nproc 16]
"""
import argparse
import os
import re
import subprocess
import tempfile

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

from multiprocessing import Pool

import numpy as np

from pp.compressibleStability import CompressibleStability, Gas

DATA = os.path.expanduser("~/Desktop/PhD/data/lst")
SEC = os.path.join(DATA, "sections")
CASES = ["d3_w13", "d1.5_w25", "d0.75_w90"]
FLOWS = ["inc", "0.05", "0.2", "0.4", "0.6", "0.8"]
RE = 1000.0
F = np.linspace(5, 400, 48)
N, YMAX, YTOP, NY = 81, 45.0, 60.0, 600
C_MIN, C_MAX, RATIO_MAX = 0.15, 0.75, 0.9
XMAX = 500.0  # legacy domains end at W + 500
XMAX_OVERRIDE = {("0.6", "d1.5_w25"): 300.0, ("0.6", "d0.75_w90"): 215.0}
LEGACY_CHK = {("0.6", "d3_w13"): "mesh_20.chk", ("0.6", "d1.5_w25"): "mesh_7.chk",
              ("0.6", "d0.75_w90"): "mesh_27.chk", ("0.8", "d3_w13"): "mesh_12.chk"}
INC_CHK = {"d3_w13": ("d3_w13_randIC", "mesh_89.chk"), "d1.5_w25": ("d1.5_w25", "mesh_30.chk"),
           "d0.75_w90": ("d0.75_w90", "mesh_23.chk")}


def dw(case):
    d, w = re.match(r"d([0-9.]+)_w([0-9.]+)", case).groups()
    return float(d), float(w)


def stations(case, flow):
    d, w = dw(case)
    xmax = XMAX_OVERRIDE.get((flow, case), XMAX)
    xg = np.linspace(0, w, int(np.ceil(w / 4)) + 2)[1:-1]  # inside the gap, edges excluded
    x = np.concatenate((np.arange(-70, 0, 5.0), xg, w + np.array([0.5, 1, 2, 3.5, 5, 7.5]),
                        np.arange(w + 10, w + 100, 5.0), np.arange(w + 100, xmax + 1e-9, 20.0)))
    return x[x <= xmax]


def source(case, flow):
    if flow == "inc":
        c, chk = INC_CHK[case]
        return f"incGapRe1000/baseflow/dns/{c}", chk
    rel = f"comGapRe1000Ma{flow}/baseflow/dns/{case}"
    if (flow, case) in LEGACY_CHK:
        return rel, LEGACY_CHK[(flow, case)]
    line = subprocess.run(["ssh", "hpc", f"grep -h '^{flow} {case} ' ~/Desktop/PhD/runs/pp_steady_cases.txt"],
                          capture_output=True, text=True, check=True).stdout.split()
    return rel, line[2]


def extract(job):
    case, flow = job
    out = os.path.join(SEC, f"{flow}_{case}.npz")
    if os.path.isfile(out):
        return f"{flow} {case}: exists"
    d, w = dw(case)
    rel, chk = source(case, flow)
    x = stations(case, flow)
    s = np.linspace(0, 1, NY) ** 2
    rows = []
    for xi in x:
        y0 = -d if 0 < xi < w else 0.0
        y = y0 + (YTOP - y0) * s
        y[0] = y0 + 1e-6  # just inside the domain
        rows.append(np.column_stack((np.full(NY, xi), y)))
    pts = np.vstack(rows)
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["rsync", "-a", f"hpc:Desktop/PhD/runs/{rel}/mesh.xml", f"hpc:Desktop/PhD/runs/{rel}/{chk}",
                        tmp + "/"], check=True)
        p = os.path.join(tmp, "q.pts")
        np.savetxt(p, pts, fmt="%.10f", comments="",
                   header='<?xml version="1.0" encoding="utf-8" ?>\n<NEKTAR>\n<POINTS DIM="2" FIELDS="">',
                   footer="</POINTS>\n</NEKTAR>")
        r = subprocess.run(["FieldConvert", "-f", "-m", f"interppoints:fromxml={tmp}/mesh.xml:fromfld={tmp}/{chk}:topts={p}",
                            os.path.join(tmp, "v.pts")], cwd=tmp, capture_output=True, text=True)
        txt = open(os.path.join(tmp, "v.pts")).read() if os.path.isfile(os.path.join(tmp, "v.pts")) else ""
        if r.returncode or not txt:
            return f"{flow} {case}: FieldConvert failed {r.stderr[-300:]}"
    fields = re.search(r'FIELDS="([^"]*)"', txt).group(1).split(",")
    body = re.search(r"<POINTS[^>]*>(.*?)</POINTS>", txt, re.S).group(1)
    v = np.array([l.split() for l in body.strip().split("\n")], float)
    assert len(v) == len(pts), (len(v), len(pts))
    data = {f.strip(): v[:, 2 + i].reshape(len(x), NY) for i, f in enumerate(fields)}
    np.savez(out, x=x, y=pts[:, 1].reshape(len(x), NY), src=f"{rel}/{chk}", **data)
    return f"{flow} {case}: {len(x)} stations from {rel}/{chk}, fields {list(data)}"


def gas_of(flow):
    return Gas(mach=0.0, viscosity="constant") if flow == "inc" else Gas(mach=float(flow))


def select(a, om):
    with np.errstate(all="ignore"):
        c = om / a.real
    ok = (a.real > 1e-6) & (c > C_MIN) & (c < C_MAX) & (np.abs(a.imag) < RATIO_MAX * a.real) & \
         (np.abs(a) < 1.5 * om / C_MIN)
    return complex(a[ok][np.argmin(a[ok].imag)]) if ok.any() else complex(np.nan, np.nan)


def solve_station(args):
    flow, case, i = args
    z = np.load(os.path.join(SEC, f"{flow}_{case}.npz"))
    y = z["y"][i]
    gas = gas_of(flow)
    try:
        if flow == "inc":
            cs = CompressibleStability(y, z["u"][i], np.ones_like(y), RE, gas, n=N, ymax=YMAX)
        else:
            cs = CompressibleStability.from_conservative(y, z["rho"][i], z["rhou"][i], z["E"][i], RE, gas,
                                                         rhov=z["rhov"][i], n=N, ymax=YMAX)
    except ValueError as e:
        return np.full(F.size, np.nan + 1j * np.nan), str(e)[:80]
    out = []
    for f in F:
        om = f * 1e-6 * RE
        a, _ = cs.spatial(om, vectors=False)
        out.append(select(a, om))
    return np.array(out), ""


def solve(nproc):
    jobs, keys = [], []
    for flow in FLOWS:
        for case in CASES:
            if os.path.isfile(os.path.join(DATA, f"gap_{flow}_{case}.npz")):
                continue
            nx = len(np.load(os.path.join(SEC, f"{flow}_{case}.npz"))["x"])
            keys.append((flow, case, nx))
            jobs += [(flow, case, i) for i in range(nx)]
    with Pool(nproc) as p:
        res = p.map(solve_station, jobs, chunksize=1)
    k = 0
    for flow, case, nx in keys:
        r = res[k:k + nx]
        k += nx
        alpha = np.array([a for a, _ in r])  # (nx, nF)
        bad = [(i, e) for i, (_, e) in enumerate(r) if e]
        x = np.load(os.path.join(SEC, f"{flow}_{case}.npz"))["x"]
        np.savez(os.path.join(DATA, f"gap_{flow}_{case}.npz"), x=x, F=F, alpha=alpha)
        print(f"{flow} {case}: {nx} stations, {len(bad)} failed {[(x[i], e) for i, e in bad][:3]}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--extract", action="store_true")
    ap.add_argument("--solve", action="store_true")
    ap.add_argument("--nproc", type=int, default=16)
    a = ap.parse_args()
    os.makedirs(SEC, exist_ok=True)
    if a.extract:
        with Pool(4) as p:
            for msg in p.imap(extract, [(c, f) for f in FLOWS for c in CASES]):
                print(msg, flush=True)
    if a.solve:
        solve(a.nproc)
