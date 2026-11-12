"""
Spatial neutral curves (alpha_i = 0, 2D waves, local parallel-flow theory) at the Mach numbers of the compressible
DNS and incompressible, with the thermodynamics of the Nektar++ sessions (gamma 1.4, Pr 0.72, Sutherland with
Tref = 225 K, isothermal wall at T_inf).

    (a) flat plate, similarity profiles, on (Re_delta*, F): scripts/lstFlatPlateMach.py -> data/lst/flatPlateMach.npz
        The band is the Re_delta* range of the DNS domain, x in [-70, 925] (Re_delta* = 1000 sqrt(1 + x C^2/1000)).
    (b-d) d3_w13, d1.5_w25, d0.75_w90 on (x, F), from the frozen DNS base-flow profiles:
        scripts/lstGapMach.py -> data/lst/gap_<flow>_<case>.npz. Thin dashed: the flat-plate curve of the same flow,
        mapped to x. Grey band: the gap. Ma0.6 (all) and Ma0.8 d3_w13 use legacy, not fully converged base flows;
        Ma0.6 d1.5_w25 / d0.75_w90 stop at x = 300 / 215 (see lstGapMach.py).

F = omega nu / U^2 x 1e6. Where no admissible mode exists (phase speed outside (0.15, 0.75)) the point counts as
stable.

Output: images/neutralCurvesMach.pdf
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from pp.colors import colors
from pp.figureFrame import CM_TO_INCH, jfm_ticks

plt.style.use("plots/style/jfm.mplstyle")

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../..")
DATA = os.path.join(ROOT, "data/lst")
FIG = os.path.join(ROOT, "images/neutralCurvesMach.pdf")
K = 1.7207876573  # Blasius delta* / sqrt(nu x / U)
RE_LE = 1000.0
FLOWS = ["inc", "0.05", "0.2", "0.4", "0.6", "0.8"]
CASES = ["d3_w13", "d1.5_w25", "d0.75_w90"]
_cm = plt.get_cmap("plasma")
COL = {"inc": "black", **{M: _cm(v) for M, v in zip(FLOWS[1:], np.linspace(0.1, 0.85, 5))}}
LAB = {"inc": "incompressible", **{M: rf"$\mathrm{{Ma}} = {M}$" for M in FLOWS[1:]}}


def x_of_re(re):
    return (np.asarray(re) ** 2 / RE_LE**2 - 1.0) * RE_LE / K**2


def re_of_x(x):
    return RE_LE * np.sqrt(1.0 + np.asarray(x) * K**2 / RE_LE)


def neutral(ax, X, Y, alpha, **kw):
    ai = np.where(np.isfinite(alpha.imag), alpha.imag, 1.0)  # no admissible mode -> stable
    return ax.contour(X, Y, ai, levels=[0.0], **kw)


def main():
    fl = np.load(os.path.join(DATA, "flatPlateMach.npz"))
    key = {0.0: "inc", 0.05: "0.05", 0.2: "0.2", 0.4: "0.4", 0.6: "0.6", 0.8: "0.8"}

    pad, w, h, gx, gy, pb, pt = 1.2, 5.2, 3.7, 0.6, 1.0, 0.95, 0.5  # width 13.4 cm (JFM 13.49)
    fw, fh = 2 * pad + 2 * w + gx, pb + 2 * h + gy + pt
    fig = plt.figure(figsize=(fw * CM_TO_INCH, fh * CM_TO_INCH))

    def axes(i, j):
        return fig.add_axes([(pad + j * (w + gx)) / fw, (pb + (1 - i) * (h + gy)) / fh, w / fw, h / fh])

    # (a) flat plate
    ax = axes(0, 0)
    R, Fg = np.meshgrid(fl["re"], fl["F"], indexing="ij")
    for k, M in enumerate(fl["mach"]):
        f = key[round(float(M), 2)]
        neutral(ax, R, Fg, fl["alpha"][k], colors=[COL[f]], linewidths=1.1)
        ax.plot([], [], color=COL[f], label=LAB[f])
    ax.axvspan(re_of_x(-70), re_of_x(925), color="0.9", zorder=0)
    ax.set_xlim(250, 4000)
    ax.set_ylim(0, 320)
    ax.set_xlabel(r"$Re_{\delta^*}$")
    ax.set_ylabel(r"$F$", rotation=0, ha="right", va="center")
    ax.set_title("(a) flat plate", fontsize=9)
    jfm_ticks(ax)
    ax.legend(loc="upper right", frameon=False, fontsize=7, handlelength=1.5)

    # (b-d) gaps
    letters = "bcd"
    for n, case in enumerate(CASES):
        i, j = divmod(n + 1, 2)
        ax = axes(i, j)
        d, wg = map(float, case[1:].split("_w"))
        ax.axvspan(0, wg, color="0.9", zorder=0)
        for f in FLOWS:
            p = os.path.join(DATA, f"gap_{f}_{case}.npz")
            if not os.path.isfile(p):
                continue
            z = np.load(p)
            X, Fg = np.meshgrid(z["x"], z["F"], indexing="ij")
            neutral(ax, X, Fg, z["alpha"], colors=[COL[f]], linewidths=1.1)
            k = [key[round(float(M), 2)] for M in fl["mach"]].index(f)
            xr = x_of_re(fl["re"])
            m = (xr >= -70) & (xr <= 520)
            if m.sum() > 2:
                Xr, Fr = np.meshgrid(xr[m], fl["F"], indexing="ij")
                neutral(ax, Xr, Fr, fl["alpha"][k][m], colors=[COL[f]], linewidths=0.6, linestyles="--")
        ax.set_xlim(-70, 500)
        ax.set_ylim(0, 400)
        ax.set_title(rf"({letters[n]}) $d = {d:g}$, $w = {wg:g}$", fontsize=9)
        ax.set_xlabel(r"$x$")
        if j == 0:
            ax.set_ylabel(r"$F$", rotation=0, ha="right", va="center")
        jfm_ticks(ax)
    fig.savefig(FIG, format="pdf")
    print(colors.OKGREEN + f"Figure saved to {FIG}" + colors.ENDC)


if __name__ == "__main__":
    main()
