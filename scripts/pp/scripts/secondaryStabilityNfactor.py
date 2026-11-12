"""
Integrated secondary amplification of the 2D waves of the d = 1.5 limit cycles at
Re = 1000, across the region III / IV boundary (w = 60, 63, 64 | 65, 66, 70).

For each width the local secondary-stability problem (pp.secondaryStability, as in
scripts/secondaryStabilityScan.py) is solved for the fundamental mode (gamma = 0) at
stations x - w = 25..X_END and a few spanwise wavenumbers, and the spatial rate
sigma = Re(sigma_t) / c is integrated downstream,

    N(x) = int_{x_0}^{x} max(sigma, 0) dx',

the amplification available to a 3D disturbance seeded at x_0 (the gap). Breakdown at
seed level A requires A e^N ~ zeta_b, i.e. N >= ln(zeta_b / A); for A = machine
precision and zeta_b = 1e-2, ln(1e-2 / 1e-16) = 32.

Output: data/stability/secondaryStabilityNfactor_d1.5.npz with widths, stations,
betas, sigma (temporal, complex) [w, station, beta], alpha, c, amp [w, station].

Run from scripts/pp: OMP_NUM_THREADS=1 .venv/bin/python scripts/secondaryStabilityNfactor.py [w ...]
(missing widths are added to an existing output file)
"""

import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from secondaryStabilityScan import OMEGA, RE, load_wave  # noqa: E402

PHD = "/home/victor/Desktop/PhD"
OUT = f"{PHD}/data/stability/secondaryStabilityNfactor_d1.5.npz"
WIDTHS = (60, 63, 64, 65, 66, 70)
STATIONS = np.arange(25, 601, 25)
BETAS = np.array([0.3, 0.5, 0.7, 2 * np.pi / 6, 1.5, 2.1])

WAVES = {}


def task(args):
    w, x0 = args
    if w not in WAVES:
        WAVES[w] = load_wave(w)
    op, L = WAVES[w].operator(x0, RE)
    return np.array([op.eigenvalues(b, 0.0)[0] for b in BETAS]), L["alpha"], L["c"], L["amp"]


def main():
    widths = [int(w) for w in sys.argv[1:]] or list(WIDTHS)
    assert all(w in OMEGA for w in widths)
    old = dict(np.load(OUT)) if os.path.exists(OUT) else None
    if old is not None:
        widths = [w for w in widths if w not in old["widths"]]
    jobs = [(w, x0) for w in widths for x0 in STATIONS]
    with Pool(15) as p:
        R = p.map(task, jobs, chunksize=1)
    shp = (len(widths), len(STATIONS))
    new = dict(widths=np.array(widths), sigma=np.array([r[0] for r in R]).reshape(shp + (len(BETAS),)),
               alpha=np.array([r[1] for r in R]).reshape(shp), c=np.array([r[2] for r in R]).reshape(shp),
               amp=np.array([r[3] for r in R]).reshape(shp))
    if old is not None:
        assert np.allclose(old["stations"], STATIONS) and np.allclose(old["betas"], BETAS)
        new = {k: np.concatenate([old[k], new[k]]) for k in new}
        o = np.argsort(new["widths"])
        new = {k: v[o] for k, v in new.items()}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    np.savez(OUT, stations=STATIONS, betas=BETAS, **new)
    print(f"Saved {OUT}: widths {list(new['widths'])}")


if __name__ == "__main__":
    main()
