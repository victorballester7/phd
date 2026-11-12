"""Convergence of a rank-1 white-noise run towards the impulse-response statistics (compressible, Ma 0.05, d1.5_w25).

Runs (HPC, comGapRe1000Ma0.05/directLinearSolver/blowingSuction/d1.5_w25, mesh R1, JacobianFreeDirect, dt = 0.2):

    meshtest_R1        impulse: q(0) = b (the vortex, weight 1), no forcing; stats = (1/t) int_0^t q q^T dt,
                       written every 28 tu (k = 1..130, t = 28 k)
    meshtest_R1_noise  q(0) = 0, forcing b a xi_n (one N(0,1) sample per time level for the whole field,
                       a^2 = 1/(0.586 dt) so that sigma^2 = 1); stats = (1/T) int_3640^{3640+T} q q^T dt,
                       written every 28 tu (k = 131..845, T = 28 (k - 130))

The driver stores raw (uncentred) products of the linearised u, v. Theory (sigma^2 = 1):

    lim_T  (1/T) int q q^T dt  (noise)  =  int_0^inf q_imp q_imp^T dt  =  t * stats_imp(t -> inf)

For each file the script computes S(x) = int_0^ymax (<u u> + <v v>) dy on the points grid (A^2 of computeNx with
doLoo = False), scaled to the impulse normalisation, and writes

    src/comGapRe1000Ma0.05/directLinearSolver/blowingSuction/d1.5_w25/data/rank1NoiseVsImpulse.npz
        x, T_noise, S_noise (n_T x n_x), t_imp, S_imp (n_t x n_x)   [S_imp = t * stats, S_noise = stats]

Usage (from scripts/pp): .venv/bin/python scripts/rank1NoiseVsImpulse.py
"""

import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HOST = "hpc"
RUN = "Desktop/PhD/runs/comGapRe1000Ma0.05/directLinearSolver/blowingSuction/d1.5_w25"
CACHE = os.path.expanduser("~/compImpulseStudy/rank1Noise")
POINTS = os.path.expanduser("~/compImpulseStudy/points_n600.pts")
OUT = os.path.expanduser(
    "~/Desktop/PhD/src/comGapRe1000Ma0.05/directLinearSolver/blowingSuction/d1.5_w25/data/rank1NoiseVsImpulse.npz")
DT_STATS = 28.0
T0_NOISE = 3640.0
K0_NOISE = 130  # noise stats k = K0_NOISE + n cover T = 28 n
N_NOISE = sorted({int(round(v)) for v in np.geomspace(1, 715, 24)})
K_IMP = list(range(5, 131, 5))


def fetch(run, ks):
    dst = os.path.join(CACHE, run)
    os.makedirs(dst, exist_ok=True)
    need = [k for k in ks if not os.path.isdir(os.path.join(dst, f"mesh_{k}_stats.fld"))]
    files = ["mesh.xml"] if not os.path.isfile(os.path.join(dst, "mesh.xml")) else []
    files += [f"mesh_{k}_stats.fld" for k in need]
    if files:
        subprocess.run(["rsync", "-aL"] + [f"{HOST}:{RUN}/{run}/{f}" for f in files] + [dst + "/"], check=True)
    return dst


def interp(d, k):
    raw = os.path.join(d, f"pts_{k}.dat")
    if not os.path.isfile(raw):
        subprocess.run(["FieldConvert", "-f", "-m", f"interppoints:fromxml={d}/mesh.xml:fromfld={d}/mesh_{k}_stats.fld:"
                        f"topts={POINTS}", raw], cwd=d, check=True, capture_output=True)
    return raw


def integrated(raw):
    """x and S(x) = int_0^ymax (uu + vv) dy on the n_x x n_y points grid."""
    names = [v.strip() for v in re.search(r"Variables\s*=\s*(.*)", open(raw).readline()).group(1).split(",")]
    d = np.loadtxt(raw, skiprows=3)
    c = {n: d[:, i] for i, n in enumerate(names)}
    x = np.unique(c["x"])
    nx = len(x)
    y = c["y"].reshape(nx, -1)
    q = (c["uu"] + c["vv"]).reshape(nx, -1)
    return x, np.trapezoid(q, y, axis=1)


def main():
    dn = fetch("meshtest_R1_noise", [K0_NOISE + n for n in N_NOISE])
    di = fetch("meshtest_R1", K_IMP)
    jobs = [(dn, K0_NOISE + n) for n in N_NOISE] + [(di, k) for k in K_IMP]
    with ThreadPoolExecutor(6) as ex:
        raws = list(ex.map(lambda a: interp(*a), jobs))
    res = [integrated(r) for r in raws]
    x = res[0][0]
    S_noise = np.array([s for _, s in res[:len(N_NOISE)]])
    t_imp = DT_STATS * np.array(K_IMP, float)
    S_imp = np.array([s for _, s in res[len(N_NOISE):]]) * t_imp[:, None]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    np.savez(OUT, x=x, T_noise=DT_STATS * np.array(N_NOISE, float), S_noise=S_noise, t_imp=t_imp, S_imp=S_imp,
             t0_noise=T0_NOISE)
    print("wrote", OUT)
    rel = np.sqrt(S_noise / S_imp[-1])
    for T, r in zip(DT_STATS * np.array(N_NOISE), rel):
        print(f"T = {T:7.0f}: A_noise/A_imp min {r.min():.3f} max {r.max():.3f} median {np.median(r):.3f}")


if __name__ == "__main__":
    main()
