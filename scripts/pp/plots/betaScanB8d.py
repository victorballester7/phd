"""
Three-dimensional global modes against spanwise wavenumber for gaps on the 2D-stable side at
Re = 1000 (d4_w15, d3_w21, d1.75_w40), single-mode Arnoldi with b = 8d, i.e. beta_k = 2 pi k / b,
k = 1..31, from data/stability/betaScan_b8d.npz (scripts/betaScanB8dCollect.py).

Top: growth rate sigma; bottom: frequency omega. Abscissa beta d = 2 pi k / 8 (top axis: k), so the
admissible wavenumbers of all domains fall on the same positions. All Ritz values with residual
< 1e-3 are shown: squares = stationary (omega = 0), circles = oscillatory; large filled markers =
leading mode at each k, small faded = sub-dominant. Lines join modes tracked from k to k + 1
(continuity in sigma and omega). Star: maximum growth. Grey crosses: old full-domain d4_w15 scan
(tau = 1750, older base flow).

Writes images/betaScanB8d.pdf.   Usage: cd scripts/pp && .venv/bin/python plots/betaScanB8d.py
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from pp.figureFrame import FigureFrame

plt.style.use("plots/style/jfm.mplstyle")

DATA = os.path.expanduser("~/Desktop/PhD/data/stability/betaScan_b8d.npz")
OUT = os.path.expanduser("~/Desktop/PhD/images/betaScanB8d.pdf")
CASES = (("d4_w15", r"$(w,d)=(15,4)$", "C3"),
         ("d3_w21", r"$(w,d)=(21,3)$", "C0"),
         ("d1.75_w40", r"$(w,d)=(40,1.75)$", "C2"))
W, H_S, H_W, GAP = 8.0, 4.4, 2.4, 0.25  # cm: frame width, sigma / omega panel heights, gap
PAD_B = 0.9
FRAME = dict(frame_w=W, aspect_ratio=H_S / W, pad_l=1.3, pad_r=3.3, pad_b=PAD_B + H_W + GAP, pad_t=0.9)
W_STAT = 5e-4  # |omega| below this: stationary
D_OM, D_SIG = 3e-3, 1.5e-3  # tracking tolerances between neighbouring k


def track(m):
    """Split modes m (rows k, sigma, omega, res) into branches by continuity from k to k + 1."""
    branches = []
    for k in np.unique(m[:, 0]):
        taken = set()
        for row in sorted(m[m[:, 0] == k], key=lambda r: -r[1]):
            best, dist = None, np.inf
            for i, b in enumerate(branches):
                last = b[-1]
                if i in taken or last[0] != k - 1:
                    continue
                dw, ds = abs(row[2] - last[2]), abs(row[1] - last[1])
                if (row[2] < W_STAT) != (last[2] < W_STAT) or dw > D_OM or ds > D_SIG:
                    continue
                if dw / D_OM + ds / D_SIG < dist:
                    best, dist = i, dw / D_OM + ds / D_SIG
            if best is None:
                branches.append([row])
                taken.add(len(branches) - 1)
            else:
                branches[best].append(row)
                taken.add(best)
    return [np.array(b) for b in branches]


def main():
    d = np.load(DATA)
    frame = FigureFrame(**FRAME)
    axs = frame.ax
    axw = frame.add_axes_cm(FRAME["pad_l"], PAD_B, W, H_W)
    for case, label, col in CASES:
        if f"{case}_modes" not in d or len(d[f"{case}_modes"]) < 5:
            continue
        m = d[f"{case}_modes"]
        m = m[~((np.abs(m[:, 1]) < 1e-5) & (m[:, 2] < W_STAT))]  # spurious |mu| ~ 1 constant mode
        bd = 2 * np.pi * m[:, 0] / 8
        lead = np.zeros(len(m), bool)
        for k in np.unique(m[:, 0]):
            i = np.flatnonzero(m[:, 0] == k)
            lead[i[np.argmax(m[i, 1])]] = True
        for b in track(m):
            if len(b) >= 3:
                x = 2 * np.pi * b[:, 0] / 8
                axs.plot(x, 1e3 * b[:, 1], color=col, lw=0.7, alpha=0.7, zorder=2)
                axw.plot(x, b[:, 2], color=col, lw=0.7, alpha=0.7, zorder=2)
        for stat, mk in ((True, "s"), (False, "o")):
            for ld, kw in ((False, dict(ms=2.2, alpha=0.4)), (True, dict(ms=3.8))):
                sel = ((m[:, 2] < W_STAT) == stat) & (lead == ld)
                for ax, y in ((axs, 1e3 * m[:, 1]), (axw, m[:, 2])):
                    ax.plot(bd[sel], y[sel], ls="none", marker=mk, color=col, mec=col if not ld else "k",
                            mew=0.3, zorder=3 + ld, **kw)
        j = np.argmax(m[:, 1])
        axs.plot(bd[j], 1e3 * m[j, 1], marker="*", ms=9, mfc=col, mec="k", mew=0.5, zorder=6)
        axw.plot(bd[j], m[j, 2], marker="*", ms=9, mfc=col, mec="k", mew=0.5, zorder=6)
        axs.annotate(rf"$\lambda_z={8 / m[j, 0]:.2f}d$", (bd[j], 1e3 * m[j, 1]), xytext=(8, 3),
                     textcoords="offset points", ha="left", va="bottom", fontsize=7)
        axs.plot([], [], color=col, lw=0.7, marker="o", ms=3.8, mec="k", mew=0.3, label=label)
    if "old_d4_w15_lzd" in d and len(d["old_d4_w15_lzd"]):
        axs.plot(2 * np.pi / d["old_d4_w15_lzd"], 1e3 * d["old_d4_w15_sigma"], ls="none", marker="x",
                 ms=3.5, color="0.45", zorder=2, label=r"$(15,4)$, earlier" "\n" r"scan, old base flow")
    axs.plot([], [], ls="none", marker="s", ms=3.8, mfc="w", mec="k", mew=0.5, label=r"stationary, $\omega=0$")
    axs.plot([], [], ls="none", marker="o", ms=3.8, mfc="w", mec="k", mew=0.5, label="oscillatory")
    axs.plot([], [], ls="none", marker="o", ms=2.2, mfc="0.6", mec="0.6", alpha=0.5, label="sub-dominant")
    xmax = 2 * np.pi * 31 / 8 + 0.4
    for ax in (axs, axw):
        ax.set_xlim(0, xmax)
    axs.axhline(0, color="k", lw=0.6, zorder=1)
    axs.set_ylim(-6, 6)
    axw.set_ylim(-0.002, 0.048)
    axs.tick_params(labelbottom=False)
    axw.set_xlabel(r"$\beta d$")
    axs.set_ylabel(r"$10^3\sigma$", rotation=0, ha="right", va="center")
    axw.set_ylabel(r"$\omega$", rotation=0, ha="right", va="center")
    top = axs.secondary_xaxis("top", functions=(lambda x: 8 * x / (2 * np.pi), lambda kk: 2 * np.pi * kk / 8))
    top.set_xticks([1, 4, 8, 12, 16, 20, 24, 28, 31])
    top.set_xlabel(r"$k$ ($b=8d$)")
    axs.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=7, borderaxespad=0)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    frame.fig.savefig(OUT)
    print("->", OUT)


if __name__ == "__main__":
    main()
