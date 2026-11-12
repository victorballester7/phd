"""Integrated secondary amplification N(x) for w = 60..70 and its calibration on the DNS.

Diagnostic for docs/notes/20261108_secondaryStability.md, from
data/stability/secondaryStabilityNfactor_d1.5.npz (scripts/secondaryStabilityNfactor.py).
Run from scripts/pp: .venv/bin/python plots/secondaryNfactor.py
Writes docs/notes/img/secondaryStability/nfactor.png; prints the calibration
N_LSA(x_b) = m ln(zeta_b / A) + b0 and the seed amplitude each width needs to break down.
Panels: wave amplitude a(x); sigma(x) (max over beta, and beta = 2 pi/6); N(x); N at the
measured x_b against ln(zeta_b / A).
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(HERE, "../../../data/stability/secondaryStabilityNfactor_d1.5.npz"))
OUT = os.path.join(HERE, "../../../docs/notes/img/secondaryStability/nfactor.png")
st, b = d["stations"], d["betas"]
COL = {60: "#9ecae1", 63: "#4292c6", 64: "#08519c", 65: "#fdae6b", 66: "#e6550d", 70: "#7f2704"}
# measured breakdown x_b - w (snapshot means) against A, from xbSnapshots / the paper
XB = {65: {1e-3: 43, 1e-5: 120, 1e-9: 233, 1e-17: 482}, 66: {1e-3: 42, 1e-5: 113, 1e-9: 211, 1e-17: 402},
      70: {1e-3: 43, 1e-5: 98, 1e-9: 187, 1e-17: 378}}


def nfactor(iw, beta=None):
    s = d["sigma"][iw].real / d["c"][iw][:, None]           # spatial rate [station, beta]
    s[d["amp"][iw] > 0.5] = 0.0                               # fit failed where the wave has vanished
    s = s.max(1) if beta is None else s[:, np.argmin(np.abs(b - beta))]
    x = np.r_[0, st]
    return s, x, np.r_[0, np.cumsum(np.maximum(np.r_[s[0], 0.5 * (s[1:] + s[:-1])], 0) * np.diff(x))]


fig, ax = plt.subplots(1, 4, figsize=(19, 4.2))
Ls, Ns = [], []
for iw, w in enumerate(d["widths"]):
    c = COL[int(w)]
    bad = d["amp"][iw] > 0.5
    smax, x, N = nfactor(iw)
    s1, _, N1 = nfactor(iw, 2 * np.pi / 6)
    ax[0].plot(st[~bad], d["amp"][iw][~bad], color=c, label=f"w={w}")
    ax[1].plot(st, smax, color=c)
    ax[1].plot(st, s1, color=c, ls="--", lw=0.8)
    ax[2].plot(x, N, color=c)
    ax[2].plot(x, N1, color=c, ls="--", lw=0.8)
    for A, xb in XB.get(int(w), {}).items():
        Ls.append(np.log(1e-2 / A)); Ns.append(np.interp(xb, x, N))
        ax[3].plot(Ls[-1], Ns[-1], "o", color=c)
    if int(w) in XB:
        ax[3].plot([], [], "o", color=c, label=f"w={w}")
m, b0 = np.polyfit(Ls, Ns, 1)
print(f"calibration: N_LSA(x_b) = {m:.2f} ln(zeta_b/A) + {b0:.2f}")
for iw, w in enumerate(d["widths"]):
    smax, x, N = nfactor(iw)
    print(f"w{w}: N_max(x-w<=600) = {N[-1]:.1f} -> breakdown needs A >~ {1e-2 * np.exp(-(N[-1] - b0) / m):.1e}"
          + ("  (lower bound on N: still growing at 600)" if smax[-1] > 0.01 else ""))
for v, lab in ((np.log(1e-2 / 1e-16), "ln(1e-2/eps)"), (np.log(1e-2 / 1e-5), "ln(1e-2/1e-5)")):
    ax[2].axhline(v, color="0.5", lw=0.6, ls=":")
    ax[2].text(590, v, lab, ha="right", va="bottom", fontsize=7)
ax[3].plot([0, 40], [0, 40], "k:", lw=0.6)
ax[3].plot([0, 40], [b0, b0 + 40 * m], "k--", lw=0.6)
ax[0].set_yscale("log"); ax[0].set_ylabel("wave amplitude a = max 2|u_1|"); ax[0].legend(fontsize=7)
ax[1].axhline(0, color="k", lw=0.5); ax[1].set_ylabel("sigma (solid: max over beta; dashed: beta=2pi/6)")
ax[2].set_ylabel("N(x) = int max(sigma, 0) dx  from the gap")
ax[3].set_xlabel("ln(zeta_b / A), zeta_b = 1e-2 (A = 1e-17 for unforced)")
ax[3].set_ylabel("LSA N at measured x_b (max over beta)")
ax[3].legend(fontsize=7)
for a in ax[:3]:
    a.set_xlabel("x - w")
fig.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=110)
print("saved", OUT)
