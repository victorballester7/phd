"""
Amplitude dependence of the secondary instability (Herbert 1988 style) of the 2D waves
of the d = 1.5 limit cycles at Re = 1000, w = 60, 63, 64 | 65, 66, 70.

At each station x - w the local wave (pp.secondaryStability.LocalWave, as in
scripts/secondaryStabilityScan.py) is rescaled to a prescribed amplitude a = max_y 2|u_1|
with the weakly nonlinear scaling (power=True: harmonic n ~ s^|n|, mean-flow distortion
~ s^2, s = a / a_DNS), keeping its shape, and the leading eigenvalue is computed for the
fundamental (gamma = 0) and subharmonic (gamma = 1/2) modes over a range of beta. The
spatial rate sigma = Re(sigma_t) / c gives the neutral curves in the (beta, a) plane and
the critical amplitude a_c = min over beta of the neutral amplitude, to compare with the
actual amplitude a_DNS(x) of each wave. Herbert's A (max rms of u') is a / sqrt 2.

Stations are kept only where the DNS wave is resolved (3e-3 < a_DNS < 0.5), and
amplitudes only up to S_MAX a_DNS (the shape of a weak wave is not reliable far beyond it).
a = 0 gives the oblique primary modes of U_BF (alpha + gamma alpha, beta).

Output: data/stability/secondaryStabilityAmplitude_d1.5.npz with widths, stations, betas,
gammas, amps; sigma[w, station, gamma, amp, beta] (spatial; nan where skipped), the
same at a = a_DNS (sigma_dns[w, station, gamma, beta]), a_dns, alpha, c [w, station].

Run from scripts/pp: OMP_NUM_THREADS=1 .venv/bin/python scripts/secondaryStabilityAmplitude.py
(about 1 h on the laptop, 15 workers; ~14 s per eigenvalue problem).
"""

import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from secondaryStabilityScan import RE, load_wave  # noqa: E402

PHD = "/home/victor/Desktop/PhD"
OUT = f"{PHD}/data/stability/secondaryStabilityAmplitude_d1.5.npz"
WIDTHS = (60, 63, 64, 65, 66, 70)
STATIONS = np.array([25, 50, 100, 200, 300, 450])
BETAS = np.array([0.2, 0.3, 0.45, 0.65, 0.9, 1.25, 1.75, 2.5])
GAMMAS = (0.0, 0.5)
AMPS = np.array([0.0, 0.0025, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16])
S_MAX = 4.0
A_MIN, A_MAX = 3e-3, 0.5

WAVES = {}


def _wave(w):
    if w not in WAVES:
        WAVES[w] = load_wave(w)
    return WAVES[w]


def info(args):
    w, x0 = args
    L = _wave(w).at(x0)
    return L["amp"], L["alpha"], L["c"]


def task(args):
    w, x0, a, a_dns = args  # a = None: the DNS amplitude
    op, L = _wave(w).operator(x0, RE, scale=1.0 if a is None else a / a_dns, power=True)
    return np.array([[op.eigenvalues(b, g)[0].real for b in BETAS] for g in GAMMAS]) / L["c"]


def main():
    ws = [(w, x0) for w in WIDTHS for x0 in STATIONS]
    with Pool(15) as p:
        I = dict(zip(ws, p.map(info, ws)))
        ok = {k: A_MIN < v[0] < A_MAX for k, v in I.items()}
        jobs = [(w, x0, a, I[w, x0][0]) for (w, x0) in ws if ok[w, x0]
                for a in [None, *AMPS] if a is None or a <= S_MAX * I[w, x0][0]]
        print(f"{sum(ok.values())} stations, {len(jobs)} x {len(BETAS) * len(GAMMAS)} eigenvalue problems", flush=True)
        R = dict(zip([j[:3] for j in jobs], p.map(task, jobs, chunksize=1)))
    nw, ns, ng, na, nb = len(WIDTHS), len(STATIONS), len(GAMMAS), len(AMPS), len(BETAS)
    sig = np.full((nw, ns, ng, na, nb), np.nan)
    sig_dns = np.full((nw, ns, ng, nb), np.nan)
    for iw, w in enumerate(WIDTHS):
        for ix, x0 in enumerate(STATIONS):
            if (w, x0, None) in R:
                sig_dns[iw, ix] = R[w, x0, None]
            for ia, a in enumerate(AMPS):
                if (w, x0, a) in R:
                    sig[iw, ix, :, ia] = R[w, x0, a]
    get = lambda k: np.array([[I[w, x0][k] for x0 in STATIONS] for w in WIDTHS])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    np.savez(OUT, widths=WIDTHS, stations=STATIONS, betas=BETAS, gammas=GAMMAS, amps=AMPS, s_max=S_MAX,
             sigma=sig, sigma_dns=sig_dns, a_dns=get(0), alpha=get(1), c=get(2))
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
