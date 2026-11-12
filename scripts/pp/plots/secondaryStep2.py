"""How LocalWave gets alpha, the harmonics u_n and the amplitude a (step 2 of the secondary analysis).

Diagnostic for docs/notes/20261108_secondaryStability.md.
Run from scripts/pp: .venv/bin/python plots/secondaryStep2.py [W=65] [X0=150]
Writes docs/notes/img/secondaryStability/step2_w<W>.png.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../scripts"))
from secondaryStabilityScan import load_wave  # noqa: E402

w = int(sys.argv[1]) if len(sys.argv) > 1 else 65
X0 = float(sys.argv[2]) if len(sys.argv) > 2 else 150
OUT = os.path.join(HERE, f"../../../docs/notes/img/secondaryStability/step2_w{w}.png")
L = load_wave(w)
x, y = L.x, L.y
W = L.at(X0)
al = W["alpha"]
lam = 2 * np.pi / al

fig, ax = plt.subplots(2, 3, figsize=(17, 8.5))
# (a) v' along x at the height of max |v'|, and its envelope
j = np.argmax(np.abs(L.v).max(0) * (y < 10))
an = hilbert(L.v[:, j])
ax[0, 0].plot(x, L.v[:, j], "C0", lw=0.8, label=f"v'(x) at y = {y[j]:.2f} (height of max |v'|)")
ax[0, 0].plot(x, np.abs(an), "k", lw=1, label="envelope |analytic signal|")
ax[0, 0].axvspan(X0 - lam, X0 + lam, color="C1", alpha=0.15, label="fit window x0 +- lambda")
ax[0, 0].set_xlim(0, 450); ax[0, 0].set_xlabel("x - w"); ax[0, 0].legend(fontsize=7)
ax[0, 0].set_title("(a) the wave along x: a carrier e^{i alpha x} under a slowly varying envelope", fontsize=9)
# (b) phase
ph = np.unwrap(np.angle(an))
ax[0, 1].plot(x, ph, "C0")
k = np.argmin(np.abs(x - X0))
ax[0, 1].plot(x, ph[k] + np.sign(ph[-1] - ph[0]) * al * (x - X0), "k--", lw=0.8,
              label=f"slope alpha = {al:.3f} at x0 = {X0:g}")
ax[0, 1].set_xlim(0, 450); ax[0, 1].set_xlabel("x - w"); ax[0, 1].set_ylabel("unwrapped phase of v' (rad)")
ax[0, 1].legend(fontsize=7)
ax[0, 1].set_title("(b) phase of the analytic signal: grows by 2 pi per wavelength", fontsize=9)
# (c) local wavenumber
ax[0, 2].plot(x, L.alpha_x, "C0", lw=0.8, label="alpha(x) = |d phase / dx|")
ax[0, 2].axvspan(X0 - 10, X0 + 10, color="C1", alpha=0.3, label="median over x0 +- 10")
ax[0, 2].axhline(al, color="k", ls="--", lw=0.8)
ax[0, 2].set_xlim(0, 450); ax[0, 2].set_ylim(0, 0.6); ax[0, 2].set_xlabel("x - w")
ax[0, 2].legend(fontsize=7)
ax[0, 2].set_title(f"(c) local wavenumber alpha = {al:.3f}, lambda = {lam:.1f},\n"
                   f"c = omega/alpha = {W['c']:.3f}", fontsize=9)
# (d) the fit at one height
m = np.abs(x - X0) <= lam
X = x[m] - X0
jy = np.argmax(np.abs(W["u"][1]))
uu = L.u[m, jy]


def rec(nmax):
    return W["u"][0][jy].real + sum(2 * (W["u"][n][jy] * np.exp(1j * n * al * X)).real for n in range(1, nmax + 1))


ax[1, 0].plot(X, uu, "ko", ms=3, label=f"DNS u' at y = {y[jy]:.2f}")
ax[1, 0].plot(X, rec(0) + 0 * X, "C2", label="u_0 only (mean)")
ax[1, 0].plot(X, rec(1), "C3", label="u_0 + u_1")
ax[1, 0].plot(X, rec(3), "C0--", label="u_0 + u_1 + u_2 + u_3")
ax[1, 0].set_xlabel("x - x0"); ax[1, 0].legend(fontsize=7)
ax[1, 0].set_title(f"(d) least-squares fit over two wavelengths at one height\n"
                   f"rms misfit with 3 harmonics: {np.sqrt(np.mean((uu - rec(3)) ** 2)):.1e} "
                   f"(signal {np.sqrt(np.mean(uu ** 2)):.1e})", fontsize=9)
# (e) harmonic profiles of u
jj = y <= 15
ax[1, 1].plot(W["u"][0][jj].real, y[jj], "C2", label="u_0 (mean-flow distortion)")
for n, c in ((1, "C3"), (2, "C1"), (3, "C4")):
    ax[1, 1].plot(2 * np.abs(W["u"][n][jj]), y[jj], c, label=f"2|u_{n}|")
ax[1, 1].axvline(W["amp"], color="C3", ls=":", lw=0.8)
ax[1, 1].text(W["amp"], 14, f"  a = max 2|u_1| = {W['amp']:.3f}", fontsize=8, color="C3")
ax[1, 1].set_xlabel("velocity"); ax[1, 1].set_ylabel("y"); ax[1, 1].legend(fontsize=7)
ax[1, 1].set_title("(e) the fit repeated at every y gives profiles", fontsize=9)
# (f) what enters the operator
ax[1, 2].plot(W["ubf"][jj], y[jj], "k", label="U_BF")
ax[1, 2].plot(W["ubf"][jj] + W["u"][0][jj].real, y[jj], "k--", label="U_BF + u_0")
ax[1, 2].axvline(W["c"], color="0.5", ls=":", label=f"c = {W['c']:.3f} (critical layer where U = c)")
ax[1, 2].set_xlabel("velocity"); ax[1, 2].set_ylabel("y"); ax[1, 2].legend(fontsize=7)
ax[1, 2].set_title("(f) mean flow seen by the disturbances\n(in the frame moving at c: U_BF + u_0 - c)", fontsize=9)
fig.suptitle(f"w = {w}, x0 - w = {X0:g}", fontsize=11)
fig.tight_layout()
fig.savefig(OUT, dpi=110)
print("saved", OUT, "| j_v:", y[j], "| j_u:", y[jy], "| alpha", al, "| rms misfit n=1:",
      np.sqrt(np.mean((uu - rec(1)) ** 2)), "n=3:", np.sqrt(np.mean((uu - rec(3)) ** 2)), "signal", np.sqrt(np.mean(uu ** 2)))
for n in range(4):
    print(f"max_y |u_{n}| = {np.abs(W['u'][n]).max():.4f}")
