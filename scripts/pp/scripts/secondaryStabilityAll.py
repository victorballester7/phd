"""
Local secondary stability of the 2D limit-cycle waves of all the d = 1.5 gaps at Re = 1000,
regions III (w = 60, 63, 64) and IV (w = 65, 66, 70), at their actual amplitude.

The same problem as scripts/secondaryStabilityScan.py (pp.secondaryStability, eq:secondaryEVP
of the JFM paper; inputs from load_wave), over the whole wave: stations x - w = 25..600 every
25, 19 spanwise wavenumbers (including beta_k = 2 pi k / 6, k = 1..4, of the Lz = 6
simulations) and gamma = 0 (fundamental) and 1/2 (subharmonic). Downstream of the first
station where the wave amplitude max_y 2|u_1| falls below A_DEAD the region III waves have
died and the fit returns noise, so those stations are skipped (nan).

Output: data/stability/secondaryStabilityAll_d1.5.npz with widths, stations, betas, gammas,
lam[w, station, gamma, beta] (complex temporal eigenvalue lambda, sigma_t = Re lambda),
alpha, c, amp, fit_err_u, fit_err_v (relative rms residual of the harmonic fit of u', v'),
alpha_hilbert (the older phase-slope estimate) [w, station], alive.
alpha is the wavenumber that minimises the residual of the harmonic fit (LocalWave.alpha_at):
the median Hilbert phase slope used before 2026-11-08 is off by up to 15 % for w = 70.

Run from scripts/pp: OMP_NUM_THREADS=1 .venv/bin/python scripts/secondaryStabilityAll.py
(about 1 h on the laptop, 15 workers).
"""

import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from secondaryStabilityScan import BETAS, GAMMAS, RE, load_wave  # noqa: E402

PHD = "/home/victor/Desktop/PhD"
OUT = f"{PHD}/data/stability/secondaryStabilityAll_d1.5.npz"
WIDTHS = (60, 63, 64, 65, 66, 70)
STATIONS = np.arange(25, 601, 25)
A_DEAD = 3e-3

WAVES = {}


def _wave(w):
    if w not in WAVES:
        WAVES[w] = load_wave(w)
    return WAVES[w]


def info(args):
    w, x0 = args
    L = _wave(w).at(x0)
    return L["amp"], L["alpha"], L["c"], L["fit_err"]["u"], L["fit_err"]["v"], _wave(w).alpha_at(x0, method="hilbert")


def task(args):
    w, x0, g = args
    op, _ = _wave(w).operator(x0, RE)
    return np.array([op.eigenvalues(b, g)[0] for b in BETAS])


def main():
    ws = [(w, x0) for w in WIDTHS for x0 in STATIONS]
    with Pool(15) as p:
        I = dict(zip(ws, p.map(info, ws)))
        alive = {}
        for w in WIDTHS:
            dead = False
            for x0 in STATIONS:
                dead = dead or I[w, x0][0] < A_DEAD
                alive[w, x0] = not dead
        jobs = [(w, x0, g) for (w, x0) in ws if alive[w, x0] for g in GAMMAS]
        print(f"{sum(alive.values())} stations, {len(jobs)} x {len(BETAS)} eigenvalue problems", flush=True)
        R = dict(zip(jobs, p.map(task, jobs, chunksize=1)))
    lam = np.full((len(WIDTHS), len(STATIONS), len(GAMMAS), len(BETAS)), np.nan, complex)
    for (w, x0, g), r in R.items():
        lam[WIDTHS.index(w), list(STATIONS).index(x0), GAMMAS.index(g)] = r
    get = lambda k: np.array([[I[w, x0][k] for x0 in STATIONS] for w in WIDTHS])
    alive_a = np.array([[alive[w, x0] for x0 in STATIONS] for w in WIDTHS])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    np.savez(OUT, widths=WIDTHS, stations=STATIONS, betas=BETAS, gammas=GAMMAS, lam=lam,
             amp=get(0), alpha=get(1), c=get(2), fit_err_u=get(3), fit_err_v=get(4), alpha_hilbert=get(5),
             alive=alive_a, a_dead=A_DEAD)
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
