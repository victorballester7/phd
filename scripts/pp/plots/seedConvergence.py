"""
Convergence check of the forced (seeded) quasi-3D runs of the d = 1.5 gaps at
Re = 1000 that do not break down without forcing (w = 63, 64). Not a paper
figure: it tells whether the runs are statistically stationary, so that they
can be stopped and their results added to fig:L2normPO (L2POvsX.py).

One column per run, three rows:
    1. energy of the spanwise Fourier modes k = 1, 2, 4 against time
       (EnergyFourierModes.mdl[.gz] on the HPC), with the times of the
       processed checkpoints (ticks) and the window used below (shaded);
    2. ||zeta||(x) of each processed checkpoint (thin), and the time average
       over the first and second halves of the window and over the whole of it,
       post-processed as in L2POvsX.py (mean of int w'^2 dy over the snapshots,
       square root, gaussian filter of the logarithm over SIGMA_3D stations).
       For A = 1e-5 the paper's w = 65 run with the same seed is drawn for
       reference. Converged: the two half-window averages coincide;
    3. breakdown position x_b - w of each checkpoint (first station past
       x - w = 25 where the raw ||zeta|| >= 1e-2, as breakdownStudy/xb.awk) with
       its cumulative mean +- one standard deviation; for runs that do not break
       down, the peak of ||zeta|| instead. Converged: the cumulative mean is flat.

The profiles (breakdownStudy/profiles/<case>_chk<n>.txt) are made from the HPC
checkpoints with the commands of breakdownStudy/procChk.job and prof.awk.
"""

import glob
import gzip
import os
import re

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d

import L2POvsX as L2

HPC = "/home/victor/hosts/hpc/PhD/runs/quasi3DGapRe1000inc/baseflow/dns"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../images/seedConvergence.pdf")
DT_CHK = 4000e3 / 300 * 0.007  # time between checkpoints (IO_CheckSteps = NumSteps/300, dt = 0.007)

# run: (width, seed, restart checkpoint, its time, first checkpoint of the window)
RUNS = {
    "d1.5_w63_Lz6.0_k32_seed1e-5": (63, "1e-5", 39, 10920.0, 52),
    "d1.5_w63_Lz6.0_k32_seed1e-9": (63, "1e-9", 39, 10920.0, 60),
    "d1.5_w64_Lz6.0_k32_seed1e-5": (64, "1e-5", 59, 5506.529, 77),
}
MODES = {1: "k = 1", 2: "k = 2", 4: "k = 4"}
X_START = 25.0  # x_b is searched downstream of this x - w
THRESHOLD = 1e-2


def chk_time(run, n):
    _, _, n0, t0, _ = RUNS[run]
    return t0 + (n - n0) * DT_CHK


def load_energy(case):
    rows = []
    for f in (f"{HPC}/{case}/EnergyFourierModes.mdl.gz", f"{HPC}/{case}/EnergyFourierModes.mdl"):
        if not os.path.exists(f):
            continue
        with (gzip.open if f.endswith("gz") else open)(f, "rt") as fh:
            for line in fh:
                p = line.split()
                if len(p) == 3 and not line.startswith("#"):
                    rows.append((float(p[0]), int(p[1]), float(p[2])))
    a = np.array(rows)
    return {k: a[a[:, 1] == k][np.argsort(a[a[:, 1] == k][:, 0])] for k in MODES}


def load_snapshots(case, w, first):
    """{chk: (x - w, int |u'|^2 dy, int w'^2 dy)} of the processed checkpoints >= first."""
    out = {}
    for f in glob.glob(f"{L2.PROFILES}/{case}_chk*.txt"):
        n = int(re.search(r"chk(\d+)", f)[1])
        if n >= first:
            d = np.loadtxt(f)
            out[n] = (d[:, 0] - w, np.abs(d[:, 1]), np.abs(d[:, 2]))
    return dict(sorted(out.items()))


def averaged(snaps, chks):
    """||zeta|| averaged over chks and post-processed as in L2POvsX.load_norm."""
    ww = np.mean([snaps[n][2] for n in chks], axis=0)
    return np.exp(gaussian_filter1d(np.log(np.sqrt(ww)), L2.SIGMA_3D))


def xb(x, ww):
    k = (x >= X_START) & (np.sqrt(ww) >= THRESHOLD)
    return x[k][0] if k.any() else np.nan


def main():
    fig, axs = plt.subplots(3, len(RUNS), figsize=(5 * len(RUNS), 10.5), constrained_layout=True)
    for j, (case, (w, seed, _, _, first)) in enumerate(RUNS.items()):
        snaps = load_snapshots(case, w, first)
        chks = list(snaps)
        t = np.array([chk_time(case, n) for n in chks])
        title = rf"$w={w}$, $A_z=10^{{{seed[2:]}}}$"

        # 1. modal energies
        ax = axs[0, j]
        for k, e in load_energy(case).items():
            ax.semilogy(e[:, 0], e[:, 2], lw=0.8, label=MODES[k])
            t_end = e[-1, 0]
        ax.axvspan(chk_time(case, first), t_end, color="0.9", zorder=0)
        for tt in t:
            ax.axvline(tt, ymax=0.04, color="k", lw=0.8)
        ax.set(title=title, xlabel="$t$", ylabel="modal energy")
        ax.legend(fontsize=8, loc="lower right")

        # 2. ||zeta||(x): snapshots, half-window and window averages
        ax = axs[1, j]
        x = snaps[chks[0]][0]
        for n in chks:
            ax.semilogy(x, np.sqrt(snaps[n][2]), color="0.75", lw=0.4)
        half = len(chks) // 2
        ax.semilogy(x, averaged(snaps, chks[:half]), color="C0", ls="--", lw=1.2,
                    label=f"chk {chks[0]}-{chks[half - 1]}")
        ax.semilogy(x, averaged(snaps, chks[half:]), color="C3", ls="--", lw=1.2,
                    label=f"chk {chks[half]}-{chks[-1]}")
        ax.semilogy(x, averaged(snaps, chks), color="k", lw=1.4, label=f"all ({len(chks)} instants)")
        if seed == "1e-5":
            xr, zr, _ = L2.load_zeta(65, "1e-5")
            ax.semilogy(xr, zr, color=L2.COLORS[65], ls=":", lw=1.2, label=r"$w=65$, same $A_z$ (paper)")
        ax.axhline(THRESHOLD, color="0.5", lw=0.6, ls=(0, (1, 2)))
        ax.set(xlim=L2.XLIM, ylim=(1e-12, 1), xlabel="$x-w$", ylabel=r"$\|\zeta\|$")
        ax.legend(fontsize=8, loc="lower right")

        # 3. breakdown position (or peak ||zeta||) against time, cumulative statistics
        ax = axs[2, j]
        vals = np.array([xb(*snaps[n][::2]) for n in chks])
        if np.isnan(vals).all():
            vals = np.array([np.sqrt(snaps[n][2]).max() for n in chks])
            ax.set_yscale("log")
            label, ylabel = r"max $\|\zeta\|$", r"max$_x\|\zeta\|$"
        else:
            label, ylabel = "$x_b - w$", "$x_b - w$"
        cm = np.array([vals[:i + 1].mean() for i in range(len(vals))])
        cs = np.array([vals[:i + 1].std() for i in range(len(vals))])
        ax.plot(t, vals, "o", ms=4, color="k", label=label)
        ax.plot(t, cm, color="C3", label="cumulative mean")
        ax.fill_between(t, cm - cs, cm + cs, color="C3", alpha=0.2, lw=0, label=r"$\pm$ s.d.")
        ax.set(xlabel="$t$", ylabel=ylabel)
        ax.legend(fontsize=8)
        print(f"{case}: chk {chks[0]}-{chks[-1]} ({len(chks)}), {label}: "
              + " ".join(f"{v:.3g}" for v in vals) + f"  -> mean {vals.mean():.3g} +- {vals.std():.2g}")
    fig.savefig(OUT)
    print(f"Saved {os.path.abspath(OUT)}")


if __name__ == "__main__":
    main()
