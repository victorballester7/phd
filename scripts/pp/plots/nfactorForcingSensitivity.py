"""
Sensitivity of the n factor of the (w, d) = (15, 4) gap at Re = 1000 to the
forcing used to excite the TS waves (appendix of the JFM paper).

The forcing of the linearised runs derives from the envelope

    E(x, y) = exp(-kx (x - xc)^2 - ky (y - yc)^2)

and a temporal signal s(t) (``scripts/makeNoiseForcing.py``). Six runs share
the base flow, mesh, time step and averaging window t in [5040, 24000]:

    d4_w15           W,       kx = ky = 100   (reference of the paper)
    d4_w15_envL      W,       kx = ky = 1
    d4_w15_envXL     W,       kx = ky = 0.2
    d4_w15_envStrip  W,       (kx, ky) = (0.05, 4), xc = -90
    d4_w15_msFlat500 S0,      kx = ky = 100
    d4_w15_msPink500 S5/3,    kx = ky = 100
    d4_w15_msFlat    S0,   F <= 250, kx = ky = 100
    d4_w15_msPink    S5/3, F <= 250, kx = ky = 100

W is unit-variance gaussian white noise sampled every time step dt = 0.006;
S0 and S5/3 are 50-line multisines, w_k = 0.01 k (10 <= F <= 500), with a
flat and a w^-5/3 spectrum (``noise_multisine.dat`` in the case folder). The
truncated runs keep the first 25 lines (F <= 250) with the same phases,
rescaled to unit variance.

Four figures: three panels, included side by side as subfigures, and the
legend shared by all of them underneath:

    nfactorForcing.pdf    n(x), with the Delta n window marked
    forcingRegions.pdf    e-folding contour of E for each envelope
    forcingSpectra.pdf    power spectral density of each signal
    forcingLegend.pdf     legend

Each forcing is also applied to the flat plate (``flatPlateRe1000inc/
directLinearSolver/blowingSuction_<variant>``, same names as the gap
variants), and Delta n of each gap case is computed against the flat plate
forced in the same way. n(x) and Delta n are computed exactly as for the
Delta n map (``pp.DeltaN_computation.computeNx`` / ``computeDeltaN``). The
averaging window of each flat plate run is in its ``data/avg_window.txt``.

Every panel keeps the same padding on both sides of its frame (FigureFrame
with pad_r = pad_l), so the frame is centred in its subfigure whatever the
width of its y labels; the two right panels share the frame width, so their
frames line up.
"""

import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, NullLocator

from pp.DeltaN_computation import computeDeltaN, computeNx
from pp.figureFrame import CM_TO_INCH, FigureFrame

plt.style.use("plots/style/jfm.mplstyle")

RE = 1000
BLASIUS_C = 1.7207876573
D, W = 4, 15
DT = 0.006  # time step of the linearised runs

ROOT = "/home/victor/Desktop/PhD"
CASE_ROOT = os.path.join(ROOT, "src", f"incGapRe{RE}", "directLinearSolver", "blowingSuction")
FLAT_ROOT = os.path.join(ROOT, "src", f"flatPlateRe{RE}inc", "directLinearSolver")

# Delta n window, as in pp.DeltaN_computation.computeDeltaN
AVG_OFFSET = 80.0
AVG_LENGTH = 50.0
AVG_NMARK = 6  # stations used to report the curve shifts over the window

XLIM_X = (-80.0, 300.0)
X_TICKS = [-70, 0, 100, 200, 300]
RE_TICKS = [900, 1000, 1100, 1200, 1300]
YLIM = (0.0, 2.5)
YTICKS = [0, 1, 2]
SHADE_GREY = "0.88"
SHADE_WINDOW = "#f6dfb8"  # Delta n window, kept apart from the grey gap band
COL_WINDOW_TEXT = "#805f26"  # darker shade of SHADE_WINDOW for its label
COL_TS_BAND = "#c2daf3"  # TS band of the spectra, opaque light blue
COL_TS_TEXT = "#145c94"
# symmetric paddings (pad_r = pad_l) centre each frame in its subfigure; the
# two right panels share the frame size so that their frames line up
RIGHT_FRAME = dict(frame_w=4.2, aspect_ratio=0.5, pad_b=0.8, pad_t=0.12)
LW_GAP = 1.3
LW_FLAT = 0.6
COL_FLAT = "tab:blue"

_green = plt.get_cmap("Greens")(np.linspace(0.55, 0.95, 3))
_purple = plt.get_cmap("Purples")(np.linspace(0.6, 0.95, 2))

# folder: (legend label, colour, linestyle, signal, kx, ky, xc)
CASES = {
    "d4_w15": (r"W, $\kappa=100$ (ref.)", "black", "-", "W", 100, 100, -95),
    "d4_w15_envL": (r"W, $\kappa=1$", _green[0], (0, (4, 1.5)), "W", 1, 1, -95),
    "d4_w15_envXL": (r"W, $\kappa=0.2$", _green[1], (0, (4, 1.5)), "W", 0.2, 0.2, -95),
    "d4_w15_envStrip": (r"W, strip", _green[2], (0, (1.2, 1.2)), "W", 0.05, 4, -90),
    "d4_w15_msFlat500": (r"$\mathrm{S}_0$, $\kappa=100$", _purple[0], (0, (5, 1.2, 1, 1.2)), "S0", 100, 100, -95),
    "d4_w15_msPink500": (r"$\mathrm{S}_{5/3}$, $\kappa=100$", _purple[1], (0, (5, 1.2, 1, 1.2)), "S53", 100, 100, -95),
    "d4_w15_msFlat": (r"$\mathrm{S}_0$, $F\leq250$", _purple[0], (0, (1, 1)), "S0", 100, 100, -95),
    "d4_w15_msPink": (r"$\mathrm{S}_{5/3}$, $F\leq250$", _purple[1], (0, (1, 1)), "S53", 100, 100, -95),
}
YC = 0.5


def flat_folder(name):
    """Flat plate run forced in the same way as the gap case ``name``."""
    return "blowingSuction" + name[len("d4_w15"):]


def flat_file(name):
    return os.path.join(FLAT_ROOT, flat_folder(name), "data", "pointsavg_n600.dat")


def x_to_re(x):
    return RE * np.sqrt(np.maximum(1.0 + BLASIUS_C**2 * np.asarray(x) / RE, 0.0))


def re_to_x(r):
    return ((np.asarray(r) / RE) ** 2 - 1.0) * RE / BLASIUS_C**2


def omega_to_F(omega):
    """F = 10^6 omega nu / u_inf^2, omega scaled with u_inf / delta*."""
    return 1e6 * np.asarray(omega) / RE


def delta_star(x):
    """Blasius delta* of the smooth wall, in units of delta* at the gap."""
    return np.sqrt(np.maximum(1.0 + BLASIUS_C**2 * np.asarray(x) / RE, 0.0))


def delta99(x):
    """Blasius delta_99 of the smooth wall, in units of delta* at the gap."""
    return 4.91 / BLASIUS_C * delta_star(x)


def plot_nfactor(save_path):
    frame = FigureFrame(frame_w=5.3, aspect_ratio=0.97, pad_l=0.65, pad_b=0.78, pad_t=0.72)
    fig, ax = frame.fig, frame.ax
    ax.grid(False)
    ax.axvspan(*x_to_re(np.array([0.0, W])), color=SHADE_GREY, lw=0, zorder=0)
    # every curve is the same gap, so the Delta n window is one band rather than markers on each curve
    x_win = x_to_re(np.array([W + AVG_OFFSET, W + AVG_OFFSET + AVG_LENGTH]))
    ax.axvspan(*x_win, color=SHADE_WINDOW, lw=0, zorder=0)
    ax.text(np.mean(x_win), YLIM[1] - 0.05, r"$\mathcal{X}$", color=COL_WINDOW_TEXT, ha="center", va="top", zorder=5)

    # reference gap and flat plate curves, to report how far each forcing shifts them over the window
    x_ref, n_ref = computeNx(flat_file("d4_w15"), False)
    x_gref, n_gref = computeNx(os.path.join(CASE_ROOT, "d4_w15", "data", "pointsavg_n600.dat"), False)
    xm = np.linspace(W + AVG_OFFSET, W + AVG_OFFSET + AVG_LENGTH, AVG_NMARK)

    print(f"\n{'case':16s} {'x0':>6s} {'dn':>7s} {'n-nref':>7s} {'nflat-nref':>10s} {'n-nflat(300)':>12s}  flat window")
    for name, (label, col, ls, *_) in CASES.items():
        x, n = computeNx(os.path.join(CASE_ROOT, name, "data", "pointsavg_n600.dat"), False)
        x_flat, n_flat = computeNx(flat_file(name), False)
        dn = computeDeltaN(W, x, n, x_flat, n_flat)
        window = os.path.join(FLAT_ROOT, flat_folder(name), "data", "avg_window.txt")
        window = open(window).read().strip() if os.path.isfile(window) else "reference run"
        print(f"{name:16s} {x[0]:6.1f} {dn:7.3f} {np.mean(np.interp(xm, x, n) - np.interp(xm, x_gref, n_gref)):7.3f}"
              f" {np.mean(np.interp(xm, x_flat, n_flat) - np.interp(xm, x_ref, n_ref)):10.3f}"
              f" {np.interp(300, x, n) - np.interp(300, x_flat, n_flat):12.3f}  {window}")
        z = 3 if name == "d4_w15" else 2
        ax.plot(x_to_re(x_flat), n_flat, color=col, ls=ls, lw=LW_FLAT, zorder=z)
        ax.plot(x_to_re(x), n, color=col, ls=ls, lw=LW_GAP, label=label, zorder=z)

    ax.set_xlim(*x_to_re(np.array(XLIM_X)))
    ax.set_ylim(*YLIM)
    ax.set_xticks(RE_TICKS)
    ax.set_yticks(YTICKS)
    ax.set_xlabel(r"$\mbox{\textit{Re}}_\mathrm{loc}$", labelpad=2)
    ax.set_ylabel(r"$n$", rotation=0, labelpad=6)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="both", top=False, right=True)
    ax.tick_params(axis="y", which="minor", length=2)
    frame.axis_on_top()
    ax2 = ax.secondary_xaxis("top", functions=(re_to_x, x_to_re))
    ax2.set_xlabel(r"$x$", labelpad=3)
    ax2.set_xticks(X_TICKS)
    ax2.xaxis.set_minor_locator(NullLocator())

    fig.savefig(save_path, format="pdf")
    plt.close(fig)
    print(f"Saved {save_path}")


def plot_regions(save_path):
    frame = FigureFrame(pad_l=0.75, **RIGHT_FRAME)
    fig, ax = frame.fig, frame.ax
    ax.grid(False)
    xlim, ylim = (-104.0, -83.0), (0.0, 3.6)

    # boundary layer and inlet
    xs = np.linspace(-100, xlim[1], 200)
    ax.plot(xs, delta99(xs), color="0.45", lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.text(-87.0, delta99(-87.0) + 0.12, r"$\delta_{99}$", color="0.35", ha="center", va="bottom")
    ax.plot(xs, delta_star(xs), color="0.45", lw=0.8, ls=(0, (1, 1.5)), zorder=1)
    ax.text(-84.5, delta_star(-84.5) + 0.08, r"$\delta^*$", color="0.35", ha="center", va="bottom")
    ax.axvspan(xlim[0], -100.0, color="0.75", lw=0, zorder=0)
    ax.text(-102.0, 3.1, "inlet", ha="center", va="center", color="0.25")

    # e-folding contour E = 1/e: an ellipse of semi-axes kx^-1/2, ky^-1/2
    t = np.linspace(0, 2 * np.pi, 400)
    labels = {
        "d4_w15_envL": (r"$\kappa=1$", (-95.0, 1.85), "center"),
        "d4_w15_envXL": (r"$\kappa=0.2$", (-95.0, 2.95), "center"),
        "d4_w15_envStrip": (r"strip", (-90.0, 1.3), "center"),
    }
    for name, (_, col, _, sig, kx, ky, xc) in CASES.items():
        if sig != "W":
            continue
        ex, ey = xc + np.cos(t) / np.sqrt(kx), YC + np.sin(t) / np.sqrt(ky)
        ax.fill(ex, ey, color=col, alpha=0.18, lw=0, zorder=2)
        # the kappa = 100 contour is under a millimetre across, so a thinner outline keeps it hollow
        ax.plot(ex, ey, color=col, lw=0.4 if name == "d4_w15" else 1.0, ls="-", zorder=3)
        if name not in labels:
            continue
        txt, pos, ha = labels[name]
        ax.text(*pos, txt, color=col, ha=ha, va="center", zorder=5)

    # kappa = 100 is a tenth of a unit across, so it is labelled with a leader
    ax.annotate(r"$\kappa=100$", xy=(-95.0, YC + 0.1), xytext=(-97.3, 1.95), ha="right", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="black", shrinkA=1, shrinkB=1.5), zorder=6)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks([-100, -95, -90, -85])
    ax.set_yticks([0, 1, 2, 3])
    ax.set_xlabel(r"$x$", labelpad=2)
    ax.set_ylabel(r"$y$", rotation=0, labelpad=8)
    frame.jfm_ticks()
    frame.axis_on_top()
    fig.savefig(save_path, format="pdf")
    plt.close(fig)
    print(f"Saved {save_path}")


def plot_spectra(save_path):
    frame = FigureFrame(pad_l=1.25, **RIGHT_FRAME)
    fig, ax = frame.fig, frame.ax
    ax.grid(False)
    ax.set_xscale("log")
    ax.set_yscale("log")
    xlim, ylim = (3.0, 2e6), (5e-7, 0.2)

    # band of TS waves unstable on the flat plate at the gap (Re_loc = 1000), from the
    # neutral points of data/localNfactor/modes_Re1000_flat_*.npz at x = 0
    ax.axvspan(55, 130, color=COL_TS_BAND, lw=0, zorder=0)
    ax.text(np.sqrt(55 * 130), 1.5e-5, "TS", color=COL_TS_TEXT, ha="center", va="center")

    # one-sided PSD in F, normalised to unit variance: white noise is flat up to the
    # Nyquist frequency, and each multisine line of variance a_k^2/2 stands for the
    # band dF = 10 (d omega = 0.01) around it
    dF = omega_to_F(0.01)
    omega_n = np.pi / DT  # Nyquist frequency of a signal sampled every time step
    F_n = omega_to_F(omega_n)
    level_w = 1.0 / F_n
    ax.plot([xlim[0], F_n, F_n], [level_w, level_w, ylim[0]], color="black", lw=1.0, zorder=3)
    ax.text(3e3, level_w * 1.25, "W", color="black", ha="center", va="bottom")

    # the truncated signals carry the same variance on half the lines and stop at F = 250.
    # S0 then sits a factor two above its full counterpart, but S5/3 has so little variance
    # above F = 250 that it moves by a few per cent, so the truncated spectra are drawn as a
    # wide pale band under the full ones rather than as a line that would hide behind them
    for name in ("d4_w15_msFlat", "d4_w15_msPink"):
        omega, amp, _ = np.loadtxt(os.path.join(CASE_ROOT, name, "noise_multisine.dat"), unpack=True)
        ax.plot(omega_to_F(omega), amp**2 / 2 / dF, color=CASES[name][1], lw=4.0, alpha=0.35,
                solid_capstyle="butt", zorder=1)
    for name, marker, mfc in (("d4_w15_msFlat500", "o", "white"), ("d4_w15_msPink500", "s", None)):
        col = CASES[name][1]
        omega, amp, _ = np.loadtxt(os.path.join(CASE_ROOT, name, "noise_multisine.dat"), unpack=True)
        psd = amp**2 / 2 / dF  # integrates to one
        F = omega_to_F(omega)
        ax.plot(F, psd, color=col, lw=0.8, zorder=2)
        ax.plot(F, psd, ls="", marker=marker, markersize=1.4, color=col,
                mfc=col if mfc is None else mfc, mew=0.6, zorder=4)
    ax.text(650, 0.004, r"$\mathrm{S}_0$", color=CASES["d4_w15_msFlat500"][1], ha="left", va="center")
    ax.text(650, 1.0e-4, r"$\mathrm{S}_{5/3}$", color=CASES["d4_w15_msPink500"][1], ha="left", va="center")

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks([1e1, 1e2, 1e3, 1e4, 1e5, 1e6])
    ax.set_yticks([1e-6, 1e-4, 1e-2])
    ax.set_xlabel(r"$F$", labelpad=2)
    ax.set_ylabel(r"$\Phi$", rotation=0, labelpad=6)
    ax.tick_params(which="both", top=True, right=True)
    ax.tick_params(which="minor", length=2)
    frame.axis_on_top()
    fig.savefig(save_path, format="pdf")
    plt.close(fig)
    print(f"Saved {save_path}")


def plot_legend(save_path):
    handles = [Line2D([], [], color=col, ls=ls, lw=LW_GAP, label=label)
               for label, col, ls, *_ in CASES.values()]
    handles += [Line2D([], [], color="0.3", lw=LW_GAP, label="Gap"),
                Line2D([], [], color="0.3", lw=LW_FLAT, label="Flat plate")]
    fig = plt.figure(figsize=(13.5 * CM_TO_INCH, 0.85 * CM_TO_INCH))
    fig.legend(handles=handles, loc="center", ncol=5, columnspacing=0.6, handlelength=2.0, handletextpad=0.5, frameon=False)
    fig.savefig(save_path, format="pdf")
    plt.close(fig)
    print(f"Saved {save_path}")


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    out_dirs = [
        os.path.join(script_path, "../../../images"),
        os.path.join(script_path, "../../../latex/papers/jfm_incNS/Images"),
    ]
    names = ["nfactorForcing.pdf", "forcingRegions.pdf", "forcingSpectra.pdf", "forcingLegend.pdf"]
    plot_nfactor(os.path.join(out_dirs[0], names[0]))
    plot_regions(os.path.join(out_dirs[0], names[1]))
    plot_spectra(os.path.join(out_dirs[0], names[2]))
    plot_legend(os.path.join(out_dirs[0], names[3]))
    for name in names:
        shutil.copy(os.path.join(out_dirs[0], name), os.path.join(out_dirs[1], name))


if __name__ == "__main__":
    main()
