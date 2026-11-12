"""
Streamwise evolution of the self-sustained oscillation of the d = 1.5 gaps at
Re = 1000 in region III (fig:L2normPO of the JFM paper).

Figures, written to images/ and to the paper's Images/:

    L2POvsX.pdf        (a) ||eta||(x) of the unforced quasi-3D runs,
                       w = 60 ... 70, with their breakdown position (w >= 65)
    L2POvsXzeta.pdf    (b) ||zeta||(x), the spanwise velocity, of the quasi-3D
                       runs at the last saved instant, for each seed amplitude A
    L2POvsXxb.pdf      (c) breakdown position x_b - w against A, statistics
                       over the saved instants of each run
    L2POvsXLegend.pdf  legend shared by (a)-(c)

Colours are the gap width: blues for the widths whose quasi-3D flow returns
to the 2D limit cycle, oranges for those that break down. Line style is the seed amplitude A, as in the sketch
(Images/L2POsketch.tex).

(a) ||eta|| = (<int |u'|^2 dy>_t)^1/2 of the unforced quasi-3D runs that break
down (w >= 65), from the same per-snapshot profiles as (b) over the stationary
window WINDOWS; its logarithm is smoothed as ||zeta||. Thin lighter lines
continue the 2D limit cycle of the same width downstream of the point where
the quasi-3D curve leaves it (ratio off by more than DEPART).
For w <= 64 (ETA_2D) the 2D runs are drawn, smoothed linearly with SIGMA = 10
stations: their quasi-3D runs converge to the 2D limit cycle but are not
converged over the whole domain. Checked on profiles processed locally with
the commands of procChk.job:
  w63 chk 30/35/39 (t = 8400-10920): max ||zeta|| 1.3e-3, 1.3e-4, 2.0e-5, and
       ||eta||_3D / ||eta||_2D tends to 1 (0.90-1.00 at chk 39);
  w64 (the unforced run stops at t = 5507 with a 2D transient in the domain;
       continuation d1.5_w64_Lz6.0_k32_kick, a spanwise kick of the seed form
       during 200 time units, then unforced) chk 104/107/110 (t = 9710-10270):
       ||zeta|| at x - w = 100 6e-17, 2e-17, 8e-18 (round-off), and the ratio
       downstream of x - w = 600 rises 0.65 -> 0.75 -> 0.86 towards 1.
  w60: the quasi-3D run stops at t = 5507, before the limit cycle settles.
The 2D w = 60 run has not converged past x - w ~ 650 either: the norm bottoms
out at 6.5e-4 and then climbs by a factor of six to the outflow, which cannot
be boundary-layer growth since the frequency of this case lies above the
upper branch everywhere downstream. Past X_EXT the curve is continued with the
decay of the neighbouring widths w = 63 and 64: the mean of their
d(ln ||eta||)/dx, scaled to the slope of w = 60 just upstream of X_EXT so that
the curve has no kink.

(b) ||zeta|| = (<int w'^2 dy>_t)^1/2 over the window WINDOWS of each run of
SEED_RUNS, from the per-snapshot wall-normal integrals in PROFILES (z-means of
u' = u - U_2D on the pointsPO_n600 points, made with the commands of
breakdownStudy/procChk.job, integrated with prof.awk). ||zeta|| spans 15
decades, so its logarithm is smoothed (SIGMA_3D stations): smoothing it
linearly would carry the large values downstream of breakdown upstream of the
forcing.

(c) x_b - w of every saved instant (breakdownStudy/xbSnapshots.dat, made with
breakdownStudy/xb.awk: first station past x - w = 25 where the raw ||zeta||
>= 1e-2) inside the window WINDOWS where each run is statistically stationary.
Markers are the mean, bars one standard deviation. The lines are least-squares
fits of the mean x_b - w against log10 A over the forced runs of the widths
that break down without forcing, dashed where extrapolated; their slope gives
sigma = ln 10 / |slope|. The unforced runs (A = 0) are placed at the machine
epsilon. Only the widths that break down without forcing are drawn.

Running the script prints the statistics, the fits, the equivalent seed
amplitude of the unforced runs and ||eta||_2D at x_b, quoted in the text.
"""

import glob
import os
import re
import shutil

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullLocator
from scipy.ndimage import gaussian_filter1d

from pp.DeltaN_computation import computeAmplitude
from pp.figureFrame import CM_TO_INCH, FigureFrame

plt.style.use("plots/style/jfm.mplstyle")

DNS_3D = "/home/victor/Desktop/PhD/src/quasi3DGapRe1000inc/baseflow/dns"
DNS_2D = "/home/victor/Desktop/PhD/src/incGapRe1000/baseflow/dns"
XB_FILE = f"{DNS_3D}/breakdownStudy/xbSnapshots.dat"
# per-snapshot wall-normal integrals (x, int |u'|^2 dy, int w'^2 dy; y <= 150) of the
# HPC data/avg/pointsPO_n600_chk<n>.dat, made by breakdownStudy/profAll.sh (prof.awk)
PROFILES = f"{DNS_3D}/breakdownStudy/profiles"

YMAX = 150  # upper limit of the wall-normal integral
SIGMA = 10  # gaussian smoothing of the 2D ||eta|| along x, in stations (dx = 2)
SIGMA_3D = 10  # the same for ln ||eta|| and ln ||zeta|| of the quasi-3D runs
X_FORCE = 20.0  # centre of the seeding force, x - w
XLIM = (0, 900)

WIDTHS_RETURN = (60, 63, 64)  # quasi-3D flow returns to the 2D limit cycle
WIDTHS_BREAK = (65, 66, 70)  # quasi-3D flow breaks down without forcing
WIDTHS = WIDTHS_RETURN + WIDTHS_BREAK
COLORS = {
    60: "#9ecae1", 63: "#4292c6", 64: "#08519c",
    65: "#fdae6b", 66: "#e6550d", 70: "#7f2704",
}
MARKERS = {65: "D", 66: "o", 70: "s"}

# seed amplitude (0: unforced) -> line style, the same as in the sketch
SEEDS = {"1e-3": "-.", "1e-5": "-", "1e-9": "--", "0": ":"}
SEED_LABELS = {"1e-3": r"$A_z=10^{-3}$", "1e-5": r"$A_z=10^{-5}$", "1e-9": r"$A_z=10^{-9}$", "0": r"$A_z=\epsilon_\mathrm{m}$"}

# Normalisation of A. makeSeedForcingCases.py divides the spanwise and temporal sums of
# the seed force by sqrt(2) and sqrt(15) (unit rms); the paper absorbs these factors into
# A, so a run with seed amplitude a has A = a / sqrt(30). To keep A a power of ten, the
# forced runs are shown at A_SHOWN = a / 10, a factor R = 10 / sqrt(30) ~ 1.83 below the
# simulated value, using x_b = x_* - ln(A) / sigma (eq:xbSeed of the paper): x_b is moved
# downstream by ln(R) / sigma, and ||zeta|| is divided by R where it is small (linear
# response to the force) and shifted downstream with x_b where it is large.
# RESCALE = False (paper version): eq. (5.4) keeps the 1/sqrt(2) and 1/sqrt(15) factors, so
# the runs are shown at their simulated amplitudes and nothing below changes the data. With
# RESCALE = True, set SEED_LABELS to 10^-4, 10^-6 and 10^-10.
RESCALE = False
NORM = np.sqrt(2.0 * 15.0)
A_SHOWN = {s: float(s) / (10 if RESCALE else 1) for s in SEEDS if s != "0"}
LN_R = np.log(10.0 / NORM) if RESCALE else 0.0
ZETA_BLEND = (1e-4, 1e-2)  # ||zeta|| range over which the rescaling turns into the shift
SEED_RUNS = {60: ("1e-5", "1e-9"), 63: ("1e-5", "1e-9"), 64: ("1e-5", "1e-9"), 65: tuple(SEEDS), 66: tuple(SEEDS), 70: tuple(SEEDS)}

# checkpoints over which each run is statistically stationary (x_b no longer
# drifts): the unforced runs from t ~ 24000 (w65, w70) and ~ 44000 (w66) on,
# after the initial noise has decayed; the forced runs once the front has
# moved to its new position, within ~ 400 time units of the restart (w65, run
# on typhoon from the unforced chk 191: within 1, 3 and 6 checkpoints for
# A = 1e-3, 1e-5 and 1e-9; the 1e-9 run processed up to chk 215).
# w63 and w64 (restarted from the unforced chk 39 and 59, stopped on 2026-11-11):
# modal energies flat from t ~ 12100 (w63 1e-5) and ~ 6700 (w64 1e-5, x_b settles
# from chk 77 on); w63 1e-9 has no breakdown and its downstream tail (x - w > 550)
# was still decaying slowly, by less than a factor of two (plots/seedConvergence.py).
# w60 1e-5/1e-9 and w64 1e-9 (restarted from the unforced chk 59, still running,
# checked on 2026-11-12 up to chk 82, t = 7653): no breakdown. w60 responds linearly
# (zeta of the two seeds differ by exactly 1e4), its modal energies are flat from
# t ~ 6000 but the leading edge of the response reaches the outflow only at chk 71;
# over 72-82 the two half-window averages agree within ~20% (a factor of 2 at
# x - w ~ 200, slow modulation of the k = 1 energy). w64 1e-9 peaks at ||zeta|| ~ 2e-3
# (x - w ~ 370, below the 1e-2 of x_b), k = 1 energy flat from t ~ 6700 (chk 72);
# the halves agree within ~20% up to x - w ~ 750, the outflow tail is not converged.
WINDOWS = {
    (65, "1e-9"): (198, 215), (65, "1e-5"): (195, 203), (65, "1e-3"): (193, 198),
    (65, "0"): (95, 191), (66, "0"): (188, 278), (70, "0"): (140, 203),
    (66, "1e-9"): (280, 327), (66, "1e-5"): (283, 301), (66, "1e-3"): (279, 301),
    (70, "1e-9"): (205, 250), (70, "1e-5"): (207, 225), (70, "1e-3"): (207, 225),
    (63, "1e-5"): (53, 71), (63, "1e-9"): (60, 70), (64, "1e-5"): (77, 90),
    (60, "1e-5"): (72, 82), (60, "1e-9"): (72, 82), (64, "1e-9"): (72, 82),
}

# (a): widths drawn from the 2D runs, whose quasi-3D runs return to the 2D limit cycle
ETA_2D = WIDTHS_RETURN
DEPART = 0.03  # relative departure of the quasi-3D ||eta|| from the 2D one where the continuation starts
X_EXT = {60: 450.0}
EXT_FROM = (63, 64)
EXT_FIT = 50.0  # length over which the slope of the extrapolated curve is measured

# the unforced runs (A = 0) are placed in (c) at the double-precision machine epsilon
EPS = float(np.finfo(float).eps)

# Symmetric paddings (pad_r = pad_l) centre each frame in its subfigure; the
# left one is set by the upright y label plus the tick labels. (a), (b) fill
# half of the JFM text width (32pc = 13.49 cm); (c) shares its row with the
# 7.50 cm sketch, so it gets 5.95 cm. All frames are 3.44 cm high.
FRAME_H = 3.44
FRAME = dict(frame_w=3.86, aspect_ratio=FRAME_H / 3.86, pad_l=1.42, pad_b=0.85, pad_t=0.15)
FRAME_C = dict(FRAME, frame_w=3.61, aspect_ratio=FRAME_H / 3.61, pad_l=1.17)


def XB_STYLE(w):
    """Breakdown position: mean (marker) and one standard deviation (capped bar)."""
    return dict(fmt=MARKERS[w], ms=4, mfc=COLORS[w], mec="k", mew=0.5, ecolor="0.15",
                elinewidth=0.7, capsize=2, capthick=0.7, zorder=4)


def case_name(w, seed):
    return f"d1.5_w{w}_Lz6.0_k32" + ("" if seed == "0" else f"_seed{seed}")


def norm_y(y, q):
    iy = np.argmin(np.abs(y - YMAX))
    return np.sqrt(np.trapezoid(np.abs(q[:, :iy]), y[:iy], axis=1))


def load_2d(w, extrapolate_tail=True):
    x, L2 = computeAmplitude(f"{DNS_2D}/d1.5_w{w}/data/pointsPO_n600.dat", False, field="|u|")
    x, L2 = x - w, gaussian_filter1d(L2, SIGMA)
    if extrapolate_tail and w in X_EXT:
        L2 = extrapolate(x, L2, X_EXT[w])
    return x, L2


def extrapolate(x, L2, xe):
    """Continue ln L2 past xe with the shape of the EXT_FROM curves."""
    ln = np.mean([np.log(load_2d(w)[1]) for w in EXT_FROM], axis=0)
    slope = np.gradient(ln, x)
    ie, i0 = np.searchsorted(x, xe), np.searchsorted(x, xe - EXT_FIT)
    # smoothing the truncated curve would bend its end, so the slope is fitted
    # over [xe - EXT_FIT, xe] where the data are still clean
    k = np.polyfit(x[i0:ie + 1], np.log(L2[i0:ie + 1]), 1)[0] / slope[i0:ie + 1].mean()
    out = L2.copy()
    out[ie:] = L2[ie] * np.exp(k * (ln[ie:] - ln[ie]))
    return out


def load_norm(case, w, col, window):
    """Time-averaged (<int q dy>_t)^1/2 over the checkpoints of window, with q = |u'|^2
    (col 1, ||eta||) or w'^2 (col 2, ||zeta||)."""
    lo, hi = window
    files = [f for f in glob.glob(f"{PROFILES}/{case}_chk*.txt")
             if lo <= int(re.search(r"chk(\d+)", f)[1]) <= hi]
    d = np.array([np.loadtxt(f) for f in files])
    # |.|: the z-mean interpolated onto the points can dip below zero at the wall
    q = np.sqrt(np.abs(d[:, :, col]).mean(axis=0))
    # the checkpoints are ~1.87 periods apart, nearly phase locked to the wave, so
    # the average keeps part of its spatial modulation; a one-wavelength filter of
    # the logarithm completes the phase average
    return d[0, :, 0] - w, np.exp(gaussian_filter1d(np.log(q), SIGMA_3D)), len(files)


def load_zeta(w, seed):
    """Time-averaged ||zeta|| = (<int w'^2 dy>_t)^1/2 over the stationary window."""
    return load_norm(case_name(w, seed), w, 2, WINDOWS[(w, seed)])


def load_eta(w):
    """Time-averaged ||eta|| = (<int |u'|^2 dy>_t)^1/2 of the unforced run."""
    if w in ETA_2D:
        x, L2 = load_2d(w)
        return x, L2, 1
    return load_norm(case_name(w, "0"), w, 1, WINDOWS[(w, "0")])


def load_xb():
    """{(w, seed): array of x_b - w over the stationary window}."""
    out = {}
    for line in open(XB_FILE):
        if line.startswith("#"):
            continue
        c, n, _, xb = line.split()[:4]
        m = re.fullmatch(r"d1\.5_w(\d+)_Lz6\.0_k32(?:_seed(.+))?", c)
        if m is None:
            continue  # e.g. the w64 kick run
        key = (int(m[1]), m[2] or "0")
        lo, hi = WINDOWS.get(key, (1, 0))
        if lo <= int(n) <= hi and xb != "nan":
            out.setdefault(key, []).append(float(xb))
    return {k: np.array(v) for k, v in out.items()}


def fit_xb(xb, w):
    """Least-squares x_b - w = p0 log10 A + p1 over the forced runs of width w."""
    seeds = [s for s in SEEDS if s != "0"]
    return np.polyfit(np.log10([float(s) for s in seeds]), [xb[(w, s)].mean() for s in seeds], 1)


def sigma_of(xb, w):
    """Spatial growth rate from the slope of the fit of x_b against log10 A."""
    return np.log(10) / -fit_xb(xb, w)[0]


def shown_xb(xb):
    """x_b - w of the forced runs moved to A_SHOWN; the unforced runs are unchanged."""
    if LN_R == 0.0:
        return xb
    return {(w, s): v if s == "0" else v + LN_R / sigma_of(xb, w) for (w, s), v in xb.items()}


def fit_shown(xbs, w):
    """Least-squares x_b - w = p0 log10 A + p1 over the forced runs, at A = A_SHOWN."""
    return np.polyfit(np.log10(list(A_SHOWN.values())), [xbs[(w, s)].mean() for s in A_SHOWN], 1)


def load_zeta_shown(w, seed, sigma):
    """||zeta|| of the run, moved from its simulated amplitude to A_SHOWN."""
    x, zeta, n = load_zeta(w, seed)
    if seed == "0" or LN_R == 0.0:
        return x, zeta, n
    ln = np.log(zeta)
    scaled = ln - LN_R
    shifted = np.interp(x - LN_R / sigma, x, ln, left=np.nan)
    shifted = np.where(np.isnan(shifted), scaled, shifted)
    # smoothstep in ln ||zeta||: rescaled below ZETA_BLEND[0], shifted above ZETA_BLEND[1];
    # in between the growth is exponential at rate sigma and both coincide
    lo, hi = np.log(ZETA_BLEND)
    t = np.clip((ln - lo) / (hi - lo), 0.0, 1.0)
    t = t * t * (3 - 2 * t)
    return x, np.exp((1 - t) * scaled + t * shifted), n


def style(frame, ylabel):
    ax = frame.ax
    ax.set_xlim(*XLIM)
    ax.set_xticks([0, 250, 500, 750])
    ax.set_yscale("log")
    ax.set_xlabel(r"$x - w$", labelpad=2)
    ax.set_ylabel(ylabel, rotation=0, labelpad=4, va="center", ha="right")
    ax.grid(False)
    frame.axis_on_top()
    frame.jfm_ticks()
    # no minor ticks on the log axis: one per unit inside each decade crowds it
    ax.yaxis.set_minor_locator(NullLocator())


def save(fig, path):
    fig.savefig(path, format="pdf")
    plt.close(fig)
    print(f"Saved {path}")


def plot_eta(save_path, xb):
    frame = FigureFrame(**FRAME)
    ax = frame.ax
    for w in WIDTHS:
        x, eta, n = load_eta(w)
        ax.plot(x, eta, color=COLORS[w], zorder=3)
        print(f"(a) w{w}: {n} instants")
        if w not in ETA_2D:
            # continuation of the 2D limit cycle past the departure of the quasi-3D flow
            x2, L2 = load_2d(w)
            off = np.abs(eta / np.interp(x, x2, L2) - 1) > DEPART
            x_dep = x[np.argmax(off & (x > 100))]
            k = x2 >= x_dep
            ax.plot(x2[k], L2[k], color=COLORS[w], lw=0.6, alpha=0.6, zorder=2)
            print(f"(a) w{w}: 2D continuation from x - w = {x_dep:.0f}")
        if (w, "0") in xb:
            # breakdown position: mean and one standard deviation over the window
            v = xb[(w, "0")]
            ax.errorbar(v.mean(), np.interp(v.mean(), x, eta), xerr=v.std(), **XB_STYLE(w))
    ax.set_ylim(3e-4, 1.2)
    style(frame, r"$\|\eta\|$")
    save(frame.fig, save_path)


def plot_zeta(save_path, sigma):
    frame = FigureFrame(**FRAME)
    ax = frame.ax
    for w, seeds in SEED_RUNS.items():
        for seed in seeds:
            x, zeta, n = load_zeta_shown(w, seed, sigma.get(w))
            ax.plot(x, zeta, ls=SEEDS[seed], color=COLORS[w], zorder=2)
            print(f"(b) w{w} A={seed}: {n} instants")
    ax.axvline(X_FORCE, color="0.6", lw=0.6, zorder=1)
    ax.set_ylim(1e-17, 1.0)
    ax.set_yticks([1e-15, 1e-10, 1e-5, 1e0])
    style(frame, r"$\|\zeta\|$")
    save(frame.fig, save_path)


def plot_xb(save_path, xb):
    xbs = shown_xb(xb)
    frame = FigureFrame(**FRAME_C)
    ax = frame.ax
    a_line = np.logspace(np.log10(EPS), -1.5, 100)
    for w in WIDTHS_BREAK:
        p = fit_shown(xbs, w)
        forced = a_line >= min(A_SHOWN.values()) / 3
        ax.plot(np.polyval(p, np.log10(a_line[forced])), a_line[forced], color=COLORS[w], lw=0.8, zorder=2)
        ax.plot(np.polyval(p, np.log10(a_line[~forced])), a_line[~forced], color=COLORS[w], lw=0.8,
                ls=(0, (3, 2)), zorder=2)
        x0 = xb[(w, "0")].mean()
        print(f"w{w}: {p[0]:.1f} per decade, sigma = {np.log(10) / -p[0]:.3f}, "
              f"A_eq = 10^{(x0 - p[1]) / p[0]:.1f}, fit at eps: x_b - w = {np.polyval(p, np.log10(EPS)):.0f}")
    # machine precision, where the unforced runs (A = 0) are placed
    ax.axhline(EPS, color="0.6", lw=0.6, ls=(0, (1, 1.5)), zorder=1)
    # widths nudged apart in A so that their bars do not overlap
    nudge = {65: 10**0.5, 66: 10**0.25, 70: 10**-0.25}
    # only the widths that break down without forcing: the others have a single
    # amplitude with breakdown (w63, w64 at 1e-5), shown in (b) and quoted in the text
    for (w, seed), v in sorted(xbs.items()):
        if w not in WIDTHS_BREAK:
            print(f"x_b - w  w{w} A={seed}: mean {v.mean():.0f}, std {v.std():.0f}, {len(v)} instants (not drawn)")
            continue
        a = (EPS if seed == "0" else A_SHOWN[seed]) * nudge[w]
        ax.errorbar(v.mean(), a, xerr=v.std(), **XB_STYLE(w))
        print(f"x_b - w  w{w} A={seed}: mean {v.mean():.0f}, range {v.min():.0f}-{v.max():.0f}, "
              f"std {v.std():.0f}, {len(v)} instants")
    ax.set_xlim(*XLIM)
    ax.set_xticks([0, 250, 500, 750])
    ax.set_yscale("log")
    ax.set_ylim(EPS / 30, 0.1)
    ax.yaxis.set_major_locator(FixedLocator([EPS, 1e-10, 1e-5]))
    ax.set_yticklabels([r"$\epsilon_\mathrm{m}$", r"$10^{-10}$", r"$10^{-5}$"])
    ax.set_xlabel(r"$x_\mathrm{b} - w$", labelpad=2)
    ax.set_ylabel(r"$A_z$", rotation=0, labelpad=1, va="center", ha="right")
    ax.grid(False)
    frame.axis_on_top()
    frame.jfm_ticks()
    ax.yaxis.set_minor_locator(NullLocator())
    save(frame.fig, save_path)


def plot_legend(save_path):
    """Three blocks: widths (2 x 3, widths returning to the limit cycle on top), amplitudes (2 x 2), unforced breakdown."""
    widths = [Line2D([], [], color=COLORS[w], label=rf"$w={w}$") for w in WIDTHS]
    styles = [Line2D([], [], color="0.3", ls=ls, label=SEED_LABELS[seed]) for seed, ls in SEEDS.items()]
    # an errorbar container drawn on a throwaway axes, so the legend shows the capped bar
    mark = [plt.figure().add_subplot().errorbar(0, 0, xerr=1, **dict(XB_STYLE(66), mfc="0.6"),
                                               label=r"$x_\mathrm{b}$: mean $\pm$ s.d.")]

    def rows_first(h, ncol):
        # matplotlib fills legends column by column; reorder so the list reads row by row
        nrow = -(-len(h) // ncol)
        return [h[r * ncol + c] for c in range(ncol) for r in range(nrow) if r * ncol + c < len(h)]

    # full JFM text width, so the x_b entry can sit at its right edge; the
    # anchors are in cm from the left edge of the figure
    fig_w = 13.49
    fig = plt.figure(figsize=(fig_w * CM_TO_INCH, 1.2 * CM_TO_INCH))
    kw = dict(frameon=False, columnspacing=0.8, handlelength=2.2, handletextpad=0.4, borderaxespad=0.0)
    for handles, ncol, xc, loc in ((widths, 3, 0.145, "center left"), (styles, 2, 8.44, "center"),
                                   (mark, 1, fig_w, "center right")):
        fig.add_artist(fig.legend(handles=rows_first(handles, ncol), ncol=ncol, loc=loc,
                                  bbox_to_anchor=(xc / fig_w, 0.5), **kw))
    save(fig, save_path)


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    out_dirs = [
        os.path.join(script_path, "../../../images"),
        os.path.join(script_path, "../../../latex/papers/jfm_incNS/Images"),
    ]
    names = ["L2POvsX.pdf", "L2POvsXzeta.pdf", "L2POvsXxb.pdf", "L2POvsXLegend.pdf"]
    xb = load_xb()
    plot_eta(os.path.join(out_dirs[0], names[0]), xb)
    plot_zeta(os.path.join(out_dirs[0], names[1]), {w: sigma_of(xb, w) for w in WIDTHS_BREAK})
    plot_xb(os.path.join(out_dirs[0], names[2]), xb)
    plot_legend(os.path.join(out_dirs[0], names[3]))
    for name in names:
        shutil.copy(os.path.join(out_dirs[0], name), os.path.join(out_dirs[1], name))


if __name__ == "__main__":
    main()
