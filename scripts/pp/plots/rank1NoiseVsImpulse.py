"""Rank-1 white-noise run against the impulse-response statistics (compressible, Ma 0.05, d1.5_w25).

Reads data/rank1NoiseVsImpulse.npz (scripts/rank1NoiseVsImpulse.py) and plots
    (a) log(A_noise / A_imp) against x for several averaging windows T_a (0 if the two agree);
    (b) max over x of |log(A_noise / A_imp)| against T_a, with a T_a^{-1/2} reference, and the truncation error
        of the impulse run, max over x of |log(A_imp(t_f) / A_imp(3640))|, against its final time t_f.
log A = 0.5 log S, S = int (uu + vv) dy; no free constant (sigma^2 = 1 in the noise run).

Usage (from scripts/pp): .venv/bin/python plots/rank1NoiseVsImpulse.py [--out file.pdf]
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np

plt.style.use("plots/style/jfm.mplstyle")

NPZ = os.path.expanduser(
    "~/Desktop/PhD/src/comGapRe1000Ma0.05/directLinearSolver/blowingSuction/d1.5_w25/data/rank1NoiseVsImpulse.npz")
W = 25.0
X_MAX = W + 900  # region of interest (outlet sponge starts at w + 900)

ap = argparse.ArgumentParser()
ap.add_argument("--out", default=os.path.expanduser("~/Desktop/PhD/images/rank1NoiseVsImpulse.pdf"))
a = ap.parse_args()

d = np.load(NPZ)
x, T, Sn, ti, Si = d["x"], d["T_noise"], d["S_noise"], d["t_imp"], d["S_imp"]
keep = x <= X_MAX
x, Sn, Si = x[keep], Sn[:, keep], Si[:, keep]
err = 0.5 * np.log(Sn / Si[-1])
err_imp = np.abs(0.5 * np.log(Si / Si[-1]))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.5, 2.5), constrained_layout=True)
pick = [np.argmin(np.abs(T - v)) for v in (112, 868, 4788, 20020)]
cols = plt.cm.viridis(np.linspace(0.15, 0.85, len(pick)))
for i, c in zip(pick, cols):
    ax1.plot(x, err[i], color=c, lw=1, label=rf"$T_\mathrm{{a}}={T[i]:.0f}$")
ax1.axhline(0, color="k", lw=0.8)
ax1.axvspan(0, W, color="0.9", lw=0)
ax1.set_xlabel(r"$x$")
ax1.set_ylabel(r"$\log(\hat A_{T_\mathrm{a}}/A)$")
ax1.legend(frameon=False, fontsize=7)
ax1.set_title("(a)", loc="left")

ax2.loglog(T, np.abs(err).max(1), "o-", ms=3, lw=1, color="C0", label="noise, max over $x$")
ax2.loglog(T, np.median(np.abs(err), 1), "s-", ms=3, lw=1, color="C1", label="noise, median over $x$")
ref = np.abs(err).max(1)[0] * np.sqrt(T[0] / T)
ax2.loglog(T, ref, "k--", lw=0.8, label=r"$\propto T_\mathrm{a}^{-1/2}$")
pos = err_imp[:-1].max(1) > 0
ax2.loglog(ti[:-1][pos], err_imp[:-1].max(1)[pos], "^-", ms=3, lw=1, color="C2", label=r"impulse, against $t_\mathrm{f}$")
ax2.set_xlabel(r"$T_\mathrm{a}$, $t_\mathrm{f}$")
ax2.set_ylabel(r"error in $\log A$")
ax2.legend(frameon=False, fontsize=7)
ax2.set_title("(b)", loc="left")

fig.savefig(a.out)
print("wrote", a.out)
