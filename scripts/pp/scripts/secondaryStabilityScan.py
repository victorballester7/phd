"""
Local secondary-stability analysis of the 2D waves of the d = 1.5 limit cycles at
Re = 1000 (w = 65, 66, 70), for sec:3Ddynamics of the JFM paper (fig:secStab).

Inputs, on the pointsPO_n600 points (x stations every 2, 600 points in y <= 150):
  2D wave:   src/incGapRe1000/baseflow/dns/d1.5_w<w>/data/pointsPO_n600.dat  (u - U_BF)
  base flow: src/quasi3DGapRe1000inc/baseflow/dns/d1.5_w<w>_Lz6.0_k32/data/bf_pointsPO_n600.dat
             (baseflow2DwithW.fld interpolated on the HPC with breakdownStudy/bfPts.job)
  omega:     frequency of v at the y = 0.5 history points of the quasi-3D runs (x - w = 25-250);
             w = 60, 63, 64 from the 2D runs (all stations agree)

Output: data/stability/secondaryStability_d1.5.npz with
  scan:   sigma (complex, temporal) of the leading mode for every width, station
          x - w in STATIONS, gamma in GAMMAS and beta in BETAS, plus alpha, c and the
          peak amplitude of the fundamental of the wave, max_y 2|u_1|
  sweep:  sigma at x - w = X_SWEEP, beta = 2 pi / 6, gamma = 0, with the wave (not
          the mean flow) scaled by SCALES

Run from scripts/pp:  OMP_NUM_THREADS=1 .venv/bin/python scripts/secondaryStabilityScan.py
(about 85 min on the laptop, 15 workers).
"""

import os
from multiprocessing import Pool

import numpy as np

from pp.secondaryStability import LocalWave

PHD = "/home/victor/Desktop/PhD"
OUT = f"{PHD}/data/stability/secondaryStability_d1.5.npz"
RE = 1000.0
WIDTHS = (65, 66, 70)
OMEGA = {60: 0.1433, 63: 0.1334, 64: 0.1308, 65: 0.1282, 66: 0.1258, 70: 0.1165}
STATIONS = np.arange(50, 251, 25)
BETA_BOX = 2 * np.pi / 6.0  # fundamental spanwise wavenumber of the Lz = 6 runs
BETAS = np.unique(np.r_[np.geomspace(0.1, 4.5, 15), BETA_BOX * np.arange(1, 5)])
GAMMAS = (0.0, 0.5)
X_SWEEP = 150
SCALES = (0.1, 0.17, 0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0)


def _read(f, cols):
    d = np.loadtxt(f, skiprows=3, usecols=cols)
    x = np.unique(d[:, 0])
    ny = len(d) // len(x)
    return x, d[:ny, 1], [d[:, c].reshape(len(x), ny) for c in range(2, d.shape[1])]


def load_wave(w):
    x, y, (u, v) = _read(f"{PHD}/src/incGapRe1000/baseflow/dns/d1.5_w{w}/data/pointsPO_n600.dat", (0, 1, 2, 3))
    xb, yb, (ub,) = _read(f"{PHD}/src/quasi3DGapRe1000inc/baseflow/dns/d1.5_w{w}_Lz6.0_k32/data/bf_pointsPO_n600.dat",
                          (0, 1, 2))
    assert np.allclose(x, xb) and np.allclose(y, yb)
    return LocalWave(x - w, y, u, v, ub, OMEGA[w])


WAVES = {}


def _wave(w):
    if w not in WAVES:
        WAVES[w] = load_wave(w)
    return WAVES[w]


def scan_task(args):
    w, x0 = args
    op, L = _wave(w).operator(x0, RE)
    sig = np.array([[op.eigenvalues(b, g)[0] for b in BETAS] for g in GAMMAS])
    return sig, L["alpha"], L["c"], L["amp"]


def sweep_task(args):
    w, s = args
    op, L = _wave(w).operator(X_SWEEP, RE, scale=s)
    return op.eigenvalues(BETA_BOX, 0.0)[0], s * L["amp"]


def main():
    scan = [(w, x0) for w in WIDTHS for x0 in STATIONS]
    sweep = [(w, s) for w in WIDTHS for s in SCALES]
    with Pool(15) as p:
        R = p.map(scan_task, scan, chunksize=1)
        S = p.map(sweep_task, sweep, chunksize=1)
    shp = (len(WIDTHS), len(STATIONS))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    np.savez(OUT, widths=WIDTHS, stations=STATIONS, betas=BETAS, gammas=GAMMAS, beta_box=BETA_BOX,
             sigma=np.array([r[0] for r in R]).reshape(shp + (len(GAMMAS), len(BETAS))),
             alpha=np.array([r[1] for r in R]).reshape(shp), c=np.array([r[2] for r in R]).reshape(shp),
             amp=np.array([r[3] for r in R]).reshape(shp),
             x_sweep=X_SWEEP, scales=SCALES,
             sweep_sigma=np.array([s[0] for s in S]).reshape(len(WIDTHS), len(SCALES)),
             sweep_amp=np.array([s[1] for s in S]).reshape(len(WIDTHS), len(SCALES)))
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
