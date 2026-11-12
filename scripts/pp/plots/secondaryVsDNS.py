"""Local secondary stability at the actual wave amplitude against the 3D DNS of the d = 1.5 gaps.

Diagnostic for docs/notes/20261108_secondaryStability.md, from
data/stability/secondaryStabilityAll_d1.5.npz (scripts/secondaryStabilityAll.py) and the
per-checkpoint ||zeta|| profiles of the quasi-3D runs (breakdownStudy/profiles, the stationary
windows of plots/L2POvsX.py).
Run from scripts/pp: .venv/bin/python plots/secondaryVsDNS.py
Writes docs/notes/img/secondaryStability/lsa_vs_dns.png and lsa_all.png.

DNS: ln||zeta||(x) of every checkpoint in the window, averaged over the checkpoints; its
local slope is the least-squares slope over x0 +- HALF (markers; bars: standard deviation of
the slopes of the individual checkpoints). Only the linear stretch is used: downstream of the
forcing (x - w >= X_LO, 25 for A = 0) and with the fit window upstream of the first station
where the mean ||zeta|| exceeds ZETA_HI.
LSA: sigma = Re(lambda) / c, fundamental (gamma = 0), max over the wavenumbers of the Lz = 6
box beta_k = 2 pi k / 6, k = 1..4 (solid), and max over all beta (dotted). Left panels: the
DNS ln||zeta|| (lines) and ln||zeta||(X_LO) + int sigma dx of the LSA (dashed).
"""
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import L2POvsX as L2  # noqa: E402

D = np.load(os.path.join(HERE, "../../../data/stability/secondaryStabilityAll_d1.5.npz"))
OUT = os.path.join(HERE, "../../../docs/notes/img/secondaryStability")
HALF = 25.0
X_LO = 50.0
ZETA_HI = 10 ** -2.5
SEEDS = ("1e-3", "1e-5", "1e-9", "0")
SCOL = {"1e-3": "C3", "1e-5": "C1", "1e-9": "C2", "0": "C0"}
COL = L2.COLORS
W, ST, B, G = list(D["widths"]), D["stations"], D["betas"], list(D["gammas"])
BK = [int(np.argmin(np.abs(B - 2 * np.pi * k / 6))) for k in range(1, 5)]


def lsa(w, gamma=0.0, box=True):
    iw, ig = W.index(w), G.index(gamma)
    s = D["lam"][iw, :, ig].real / D["c"][iw][:, None]
    s = s[:, BK] if box else s
    return np.where(D["alive"][iw], np.nanmax(s, 1) if np.isfinite(s).any() else np.nan, np.nan)


def integrate(x, s, x0, xe):
    """int_{x0}^{x} sigma dx' for x0 <= x <= xe (sigma = 0 where the wave is dead)."""
    s = np.where(np.isfinite(s), s, 0.0)
    xs = np.r_[x0, ST[ST > x0]]
    ss = np.interp(xs, ST, s)
    N = np.r_[0, np.cumsum(0.5 * (ss[1:] + ss[:-1]) * np.diff(xs))]
    m = xs <= xe
    return xs[m], N[m]


def dns(w, seed):
    lo, hi = L2.WINDOWS[(w, seed)]
    files = sorted(f for f in glob.glob(f"{L2.PROFILES}/{L2.case_name(w, seed)}_chk*.txt")
                   if lo <= int(re.search(r"chk(\d+)", f)[1]) <= hi)
    Z = np.array([0.5 * np.log(np.abs(np.loadtxt(f)[:, 2]) + 1e-300) for f in files])
    x = np.loadtxt(files[0])[:, 0] - w
    m = Z.mean(0)
    x_lo = 25.0 if seed == "0" else X_LO
    x_hi = x[np.argmax((x > x_lo) & (m > np.log(ZETA_HI)))]
    xs, s, e = [], [], []
    for x0 in np.arange(x_lo + HALF, x_hi - HALF + 1, 12.5):
        k = np.abs(x - x0) <= HALF
        xs.append(x0)
        s.append(np.polyfit(x[k], m[k], 1)[0])
        e.append(np.std([np.polyfit(x[k], z[k], 1)[0] for z in Z]))
    return x, m, x_lo, x_hi, np.array(xs), np.array(s), np.array(e), len(files)


# ---- region IV: LSA against DNS ----------------------------------------------------------
fig, axs = plt.subplots(3, 2, figsize=(15, 12))
print("run          n   x_lo  x_hi   <sigma_DNS>  <sigma_LSA box>  ratio   Dln(zeta) DNS  int sigma LSA")
for row, w in zip(axs, (65, 66, 70)):
    sb, sa, s05 = lsa(w), lsa(w, box=False), lsa(w, 0.5)
    row[1].plot(ST, sb, "k-", lw=1.6, label="LSA fundamental, max over beta_k = 2 pi k/6")
    row[1].plot(ST, sa, "k:", lw=1.0, label="LSA fundamental, max over all beta")
    row[1].plot(ST, s05, "0.5", ls="--", lw=0.8, label="LSA subharmonic, max over beta_k")
    poor = D["fit_err_u"][W.index(w)] > 0.15
    row[1].plot(ST[poor], sb[poor], "ks", mfc="w", ms=6, label="LSA: wave model poor (fit residual above 15 percent)")
    for seed in SEEDS:
        x, m, x_lo, x_hi, xs, s, e, n = dns(w, seed)
        c = SCOL[seed]
        lab = "A=0" if seed == "0" else f"A={seed}"
        row[0].plot(x, m, color=c, lw=1, label=f"DNS {lab} ({n} chk)")
        xi, N = integrate(ST, sb, x_lo, x_hi + 50)
        row[0].plot(xi, np.interp(x_lo, x, m) + N, color=c, ls="--", lw=1.4)
        row[0].axvspan(x_lo, x_hi, color=c, alpha=0.05)
        row[1].errorbar(xs, s, e, fmt="o", color=c, ms=4, capsize=2, label=f"DNS {lab}")
        sl = np.interp(xs, ST, sb)
        dN_dns = np.interp(x_hi, x, m) - np.interp(x_lo, x, m)
        dN_lsa = np.interp(x_hi, *integrate(ST, sb, x_lo, 1e9))
        print(f"w{w} A={seed:5s} {n:3d} {x_lo:5.0f} {x_hi:5.0f}   {np.mean(s):.3f}        {np.mean(sl):.3f}"
              f"        {np.mean(s) / np.mean(sl):.2f}     {dN_dns:6.1f}        {dN_lsa:6.1f}")
    row[0].axhline(np.log(1e-2), color="0.5", lw=0.6, ls=":")
    row[0].set_xlim(0, 500); row[0].set_ylim(-37, 1)
    row[0].set_ylabel(r"$\ln\|\zeta\|$"); row[0].legend(fontsize=7, loc="lower right")
    row[0].set_title(f"w={w}: DNS " r"$\ln\|\zeta\|$" " (lines) and DNS(x_lo) + int sigma_LSA dx (dashed); "
                     "shaded: linear stretch", fontsize=8.5)
    row[1].set_xlim(0, 500); row[1].set_ylim(-0.02, 0.2); row[1].axhline(0, color="0.5", lw=0.5)
    row[1].set_ylabel("spatial growth rate sigma"); row[1].legend(fontsize=7, ncol=2)
    row[1].set_title(f"w={w}: LSA sigma(x) (lines) and DNS slope of " r"$\ln\|\zeta\|$" f" over +-{HALF:g} (markers)",
                     fontsize=8.5)
for a in axs[-1]:
    a.set_xlabel("x - w")
fig.tight_layout(); fig.savefig(f"{OUT}/lsa_vs_dns.png", dpi=110); plt.close(fig)

# ---- all widths: sigma(x) and N(x) ----------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(19, 4.8))
print("\nw    N_box(600)  N_all(600)   x where wave dies")
for w in W:
    iw = W.index(w)
    c = COL[w]
    axs[0].plot(ST, D["amp"][iw] * np.where(D["alive"][iw], 1, np.nan), color=c, label=f"w={w}")
    axs[1].plot(ST, lsa(w), color=c, lw=1.6)
    axs[1].plot(ST, lsa(w, box=False), color=c, ls=":", lw=1)
    xi, N = integrate(ST, np.maximum(np.nan_to_num(lsa(w)), 0), 0, 600)
    xa, Na = integrate(ST, np.maximum(np.nan_to_num(lsa(w, box=False)), 0), 0, 600)
    axs[2].plot(xi, N, color=c, lw=1.6)
    axs[2].plot(xa, Na, color=c, ls=":", lw=1)
    dead = ST[~D["alive"][iw]]
    print(f"w{w}   {N[-1]:6.1f}     {Na[-1]:6.1f}      {dead[0] if len(dead) else '-'}")
axs[0].set_yscale("log"); axs[0].set_ylabel(r"wave amplitude $a=\max_y 2|u_1|$"); axs[0].legend(fontsize=7)
axs[1].axhline(0, color="0.5", lw=0.5)
axs[1].set_ylabel("sigma (solid: max over beta_k of the box; dotted: all beta)")
axs[2].set_ylabel("N(x) = int max(sigma, 0) dx from the gap")
for v, lab in ((np.log(1e-2 / 1e-17), "ln(1e-2/1e-17)"), (np.log(1e-2 / 1e-9), "ln(1e-2/1e-9)"),
               (np.log(1e-2 / 1e-5), "ln(1e-2/1e-5)")):
    axs[2].axhline(v, color="0.5", lw=0.6, ls=":")
    axs[2].text(590, v, lab, ha="right", va="bottom", fontsize=7)
for a in axs:
    a.set_xlabel("x - w"); a.set_xlim(0, 600)
fig.suptitle("Fundamental secondary mode of the actual 2D wave, all d = 1.5 widths (curves stop where "
             "the wave has died)", fontsize=10)
fig.tight_layout(); fig.savefig(f"{OUT}/lsa_all.png", dpi=110); plt.close(fig)
print("saved", OUT)
