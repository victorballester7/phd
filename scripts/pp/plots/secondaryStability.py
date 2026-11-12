"""
Local secondary stability of the 2D waves of the d = 1.5 limit cycles at Re = 1000
(fig:secStab of the JFM paper), from data/stability/secondaryStability_d1.5.npz
(scripts/secondaryStabilityScan.py).

Figures, written to images/ and to the paper's Images/:

    secStabBeta.pdf    (a) spatial growth rate sigma = sigma_t / c against beta at
                       x - w = 150; solid fundamental (gamma = 0), dashed subharmonic
                       (gamma = 1/2); markers at the wavenumbers beta_k = 2 pi k / b
                       of the quasi-3D domain (b = 6)
    secStabX.pdf       (b) sigma at beta_1 = 2 pi / b, gamma = 0, against x - w (lines),
                       and d ln||zeta|| / dx of the quasi-3D runs forced with A = 1e-9
                       (markers), time-averaged profiles, local least-squares slope
                       over +-25 units, at stations at least 25 units upstream of the
                       mean x_b of that run
    secStabAmp.pdf     (c) temporal rate sigma_t at x - w = 150, beta_1, gamma = 0,
                       against the peak amplitude of the fundamental of the wave,
                       max_y 2|u_1|, with the wave scaled; filled markers at the
                       amplitude of the limit cycle
    secStabLegend.pdf  legend
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import L2POvsX as L2  # noqa: E402  (colours, markers, DNS windows and profiles)

from pp.figureFrame import CM_TO_INCH, FigureFrame  # noqa: E402

DATA = "/home/victor/Desktop/PhD/data/stability/secondaryStability_d1.5.npz"
FRAME = dict(frame_w=3.2, aspect_ratio=0.9, pad_l=1.0, pad_r=0.2, pad_b=0.85, pad_t=0.15)
X_BETA = 150
FIT_HALF = 25.0


def dns_rate(w, stations):
    """d ln||zeta||/dx of the A = 1e-9 run at the stations upstream of its mean x_b."""
    import glob
    import re
    lo, hi = L2.WINDOWS[(w, "1e-9")]
    files = [f for f in glob.glob(f"{L2.PROFILES}/{L2.case_name(w, '1e-9')}_chk*.txt")
             if lo <= int(re.search(r"chk(\d+)", f)[1]) <= hi]
    lnz = []
    for f in files:
        d = np.loadtxt(f)
        x = d[:, 0] - w
        lnz.append(0.5 * np.log(np.abs(d[:, 2]) + 1e-300))
    lnz = np.mean(lnz, 0)
    xb = L2.load_xb()[(w, "1e-9")].mean()
    xs, s = [], []
    for x0 in stations:
        if x0 > xb - FIT_HALF:
            continue
        m = np.abs(x - x0) <= FIT_HALF
        xs.append(x0)
        s.append(np.polyfit(x[m], lnz[m], 1)[0])
    return np.array(xs), np.array(s)


def finish(frame, xlabel, ylabel):
    ax = frame.ax
    ax.set_xlabel(xlabel, labelpad=2)
    ax.set_ylabel(ylabel, labelpad=2)
    ax.grid(False)
    frame.axis_on_top()
    frame.jfm_ticks()


def plot_beta(d, path):
    frame = FigureFrame(**FRAME)
    ax = frame.ax
    i = list(d["stations"]).index(X_BETA)
    b = d["betas"]
    for iw, w in enumerate(d["widths"]):
        c = d["c"][iw, i]
        for ig, ls in ((0, "-"), (1, "--")):
            ax.plot(b, d["sigma"][iw, i, ig].real / c, ls=ls, color=L2.COLORS[w], lw=0.9)
        kb = np.isin(np.round(b, 6), np.round(d["beta_box"] * np.arange(1, 5), 6))
        ax.plot(b[kb], d["sigma"][iw, i, 0, kb].real / c, ls="none", marker=L2.MARKERS[w], ms=3.5,
                mfc=L2.COLORS[w], mec="k", mew=0.5)
    ax.axhline(0, color="0.6", lw=0.5, zorder=1)
    ax.set_xlim(0, 4.5)
    ax.set_ylim(-0.03, 0.13)
    finish(frame, r"$\beta$", r"$\sigma$")
    L2.save(frame.fig, path)


def plot_x(d, path):
    frame = FigureFrame(**FRAME)
    ax = frame.ax
    ib = int(np.argmin(np.abs(d["betas"] - d["beta_box"])))
    for iw, w in enumerate(d["widths"]):
        s = d["sigma"][iw, :, 0, ib].real / d["c"][iw]
        ax.plot(d["stations"], s, color=L2.COLORS[w], lw=0.9)
        xs, sd = dns_rate(w, d["stations"])
        ax.plot(xs, sd, ls="none", marker=L2.MARKERS[w], ms=3.5, mfc="none", mec=L2.COLORS[w], mew=0.8)
        print(f"w{w}: LSA sigma(beta_1) = " + " ".join(f"{x:.0f}:{v:.3f}" for x, v in zip(d["stations"], s)))
        print(f"      DNS dln(zeta)/dx = " + " ".join(f"{x:.0f}:{v:.3f}" for x, v in zip(xs, sd)))
        print(f"      wave amplitude   = " + " ".join(f"{x:.0f}:{v:.3f}" for x, v in zip(d["stations"], d["amp"][iw])))
        print(f"      c                = " + " ".join(f"{x:.0f}:{v:.3f}" for x, v in zip(d["stations"], d["c"][iw])))
        smax = d["sigma"][iw].real.max(axis=(1, 2)) / d["c"][iw]
        bmax = d["betas"][d["sigma"][iw, :, 0].real.argmax(1)]
        print(f"      max over beta, gamma = " + " ".join(f"{v:.3f}" for v in smax) + "  at beta " +
              " ".join(f"{v:.2f}" for v in bmax))
    ax.set_xlim(25, 275)
    ax.set_xticks([50, 150, 250])
    ax.set_ylim(0.0, 0.13)
    finish(frame, r"$x - w$", r"$\sigma$")
    L2.save(frame.fig, path)


def plot_amp(d, path):
    frame = FigureFrame(**FRAME)
    ax = frame.ax
    one = list(d["scales"]).index(1.0)
    for iw, w in enumerate(d["widths"]):
        a, s = d["sweep_amp"][iw], d["sweep_sigma"][iw].real
        ax.plot(a, s, color=L2.COLORS[w], lw=0.9)
        ax.plot(a[one], s[one], ls="none", marker=L2.MARKERS[w], ms=3.5, mfc=L2.COLORS[w], mec="k", mew=0.5)
        k = a > 0.05
        print(f"w{w}: sigma_t ~ a^{np.polyfit(np.log(a[k]), np.log(s[k]), 1)[0]:.2f} for a > 0.05; "
              + " ".join(f"{x:.3f}:{v:.4f}" for x, v in zip(a, s)))
    ar = np.array([0.05, 0.4])
    ax.plot(ar, 0.03 * (ar / 0.15) ** 0.5, color="0.4", lw=0.6, ls=(0, (3, 2)))
    ax.text(0.2, 0.03 * (0.2 / 0.15) ** 0.5 / 1.2, r"$\propto a^{1/2}$", ha="left", va="top")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.008, 0.5)
    ax.set_ylim(0.002, 0.12)
    finish(frame, r"$a$", r"$\sigma_t$")
    L2.save(frame.fig, path)


def plot_legend(d, path):
    h = [Line2D([], [], color=L2.COLORS[w], marker=L2.MARKERS[w], ms=3.5, mfc=L2.COLORS[w], mec="k", mew=0.5,
                label=rf"$w={w}$") for w in d["widths"]]
    h += [Line2D([], [], color="0.3", ls="-", label=r"fundamental, $\gamma=0$"),
          Line2D([], [], color="0.3", ls="--", label=r"subharmonic, $\gamma=1/2$"),
          Line2D([], [], color="0.3", ls="none", marker="o", ms=3.5, mfc="none", mec="0.3", mew=0.8,
                 label=r"DNS, $A=10^{-9}$")]
    h = [h[0], h[3], h[1], h[4], h[2], h[5]]  # legend fills columns: widths on the top row
    fig = plt.figure(figsize=(13.2 * CM_TO_INCH, 0.9 * CM_TO_INCH))
    fig.legend(handles=h, loc="center", ncol=3, frameon=False, handlelength=2.0, columnspacing=1.5)
    L2.save(fig, path)


def main():
    d = dict(np.load(DATA))
    script_path = os.path.dirname(os.path.abspath(__file__))
    out_dirs = [os.path.join(script_path, "../../../images"),
                os.path.join(script_path, "../../../latex/papers/jfm_incNS/Images")]
    names = {"secStabBeta.pdf": plot_beta, "secStabX.pdf": plot_x, "secStabAmp.pdf": plot_amp,
             "secStabLegend.pdf": plot_legend}
    import shutil
    for name, f in names.items():
        f(d, os.path.join(out_dirs[0], name))
        shutil.copy(os.path.join(out_dirs[0], name), os.path.join(out_dirs[1], name))


if __name__ == "__main__":
    main()
