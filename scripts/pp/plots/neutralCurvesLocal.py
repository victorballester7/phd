"""
Neutral curves of the gap flow from the spatial local stability analysis, one
panel per family.

Same cases, same cached scans and same colours as
``plots/nfactorCurvesLocal.py``: the spatial Orr-Sommerfeld problem solved at
every station of the frozen base flow for a band of real reduced frequencies,
``scripts/predictNfactorLocal.py`` -> ``data/localNfactor/modes_Re*_*.npz``.
Where that figure integrates alpha_i into n(x), this one draws the locus

    alpha_i(x, F) = 0

on the (Re_l, F) plane, so the two can be read against each other: the lobe a
case opens here is what the n factor of that case integrates across.

Upstream of the gap every case is the same Blasius flow, and far downstream
they have all relaxed back onto it, so the curves collapse onto one common
Tollmien-Schlichting branch at both ends of the window and differ only over
and just after the gap. The fan of lobes in between is the figure.

One caveat to read the bottom of the frame with. The scan covers F = 10 to
500, and at almost every case here the lower branch of the gap lobe sits
below that, so what closes those lobes on F = 0 is not a computed neutral
point: it is the statement that the column is still amplified at the bottom
of the scanned band, and that inside the gap the inflectional shear layer is
unstable down to arbitrarily small F with its growth rate falling smoothly to
zero with the frequency. The two shallow cases whose branch is inside the
band over most of the gap are the exception, and near the upstream edge,
where even they lose it, the branch is extrapolated rather than closed --
see EXTRAPOLATE.
"""

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, NullLocator

from pp.fileManagement import extract_depth_width
from pp.localNfactor import ModeSet
from pp.figureFrame import FigureFrame

plt.style.use("plots/style/jfm.mplstyle")

RE = 1000
BLASIUS_C = 1.7207876573

ROOT = "/home/victor/Desktop/PhD"
DATA_ROOT = os.path.join(ROOT, "data", "localNfactor")

# the same scan the n factor figure reads, so the two are the same experiment
TAG = "F10-500x28_n81_y45_X700"

COL_FLAT = "tab:blue"
COL_BFS = "tab:orange"

# legend_x is the horizontal anchor of the case legend, in axes coordinates.
# It is per family because the two legends are not the same width: the w40
# block carries seven "d = x.xx" labels against six "w = xx" and a "BFS", so
# centring both on 0.5 leaves the wider one overhanging on the right. Nudging
# it left recovers that. Around 0.44 is as far as it can go before it reads as
# misaligned rather than as fitted.

FAMILIES = {
    "d1.5": dict(
        cases={
            "d1.5_w10": r"$w=10$",
            "d1.5_w15": r"$w=15$",
            "d1.5_w20": r"$w=20$",
            "d1.5_w30": r"$w=30$",
            "d1.5_w40": r"$w=40$",
            "d1.5_w50": r"$w=50$",
        },
        cmap="Greens",
        extras=[("bfs_d1.5", "BFS", COL_BFS)],
        out="neutralCurvesLocald1.5.pdf",
        legend_ncol=3,
        legend_x=0.5,
    ),
    "w40": dict(
        cases={
            "d0.25_w40": r"$d=0.25$",
            "d0.5_w40": r"$d=0.50$",
            "d0.75_w40": r"$d=0.75$",
            "d1_w40": r"$d=1.00$",
            "d1.25_w40": r"$d=1.25$",
            "d1.5_w40": r"$d=1.50$",
            "d1.75_w40": r"$d=1.75$",
        },
        cmap="Reds",
        extras=[],
        out="neutralCurvesLocalw40.pdf",
        legend_ncol=3,
        legend_x=0.5,
    ),
}

# The frame is the scanned band itself: F from 10 to 500 on the left. The
# streamwise window is shorter than the one plots/nfactorCurvesLocal.py uses,
# because the two figures need different things from it: the n factor keeps
# growing out to x = 600, whereas every lobe here is shut by x ~ 60 and the
# curves have merged back onto the flat plate branch well before x = 250.
FLIM = (0.0, 500.0)
FTICKS = [0, 150, 300, 450]
XLIM_X = (-80.0, 270.0)
# The two horizontal ladders are kept apart: Re_l ticks on the bottom edge
# only, x ticks on the top edge only, major ticks alone on both. The Re_l
# ladder starts at Re and the x one at 0, so the first tick of each marks the
# same station, the upstream gap edge. Only the ticks inside the window are
# drawn.
RE_TICKS = [900, 1000, 1100, 1200, 1300]
X_TICKS = [-70, 0, 40, 100,  200]

FIG_W_CM = 5.75
FIG_ASPECT = 0.75

LW_CASE = 0.9
LW_FLAT = 1.2

# Phase speed at which the discretised continuous spectrum piles up, the
# C_MAX of pp.localNfactor.select_mode. When no genuine discrete mode is
# admissible the selection falls back on that pile up, so c_r >= this value
# means "no instability mode was found here", not "the mode runs this fast".
C_FALLBACK = 0.70

# How far the lower branch of the lobe is repaired where the scan cannot see
# it. The scan starts at F = 10, and inside the gap the shear layer is
# unstable below that at every case here, so the raw contour leaves the lobe
# hanging open on the bottom of the band.
#
#   close_at_f0     close the lobe on F = 0 wherever the column is still
#                   amplified at the bottom of the band. That is the honest
#                   reading of what the scan found: the growth rate falls
#                   smoothly to zero with the frequency and stays positive all
#                   the way down, the inflectional shear layer being unstable
#                   at every frequency it carries.
#   EXTRAPOLATE     cases whose lower branch *is* inside the band over most of
#                   the gap but drops out of it near the upstream edge. There
#                   the branch is continued by a parabola fitted to the
#                   stations that do carry it -- a reading of the trend, not a
#                   computed neutral point. Only the two shallow cases qualify;
#                   at every other case the branch is out of band across the
#                   whole gap, so there is nothing to fit, and the step is
#                   excluded too because its branch only re-enters the band
#                   60 units downstream, far too long a run to extrapolate
#                   back over.
CLOSE_AT_F0 = True
EXTRAPOLATE = {"d0.25_w40", "d0.5_w40"}
LOW_BRANCH_FIT = 20  # stations used by that fit

# Light smoothing along x, applied to the growth rate before contouring, to
# take the staircase off the places where the contour hands over from one
# branch to another -- the shear layer mode to the Tollmien-Schlichting mode
# downstream of the gap, most visibly on the step near Re_l = 1100. It is a
# [1 2 1]/4 kernel applied SMOOTH_PASSES times *within* each of the three
# streamwise segments (before the gap, inside it, after it) and never across
# their boundaries, so the genuinely sharp edges at x = 0 and x = w survive
# untouched. Set to 0 to see the raw contour.
SMOOTH_PASSES = 2

# gap shading, drawn behind everything
SHADE_ALPHA = 0.3
SHADE_GREY = "0.88"


def gap_width(name):
    """Streamwise extent of the gap, inf for a step and 0 for the flat plate."""
    if name == "flat":
        return 0.0
    if name.startswith("bfs_"):
        return np.inf  # a step has no downstream edge
    _, w = extract_depth_width(name)
    return float(w)


def shade_gaps(ax, entries):
    """
    Mark where the gap is.

    When every case in the family has the same gap there is one band and it is
    grey, as a plain background. When the widths differ, each case colours the
    stretch its own gap adds to the one before it -- the widths are nested, all
    starting at x = 0 -- so the shading is a staircase whose steps are the
    downstream edges, and each colour appears once, at its own edge, with no
    translucent bands stacked on top of each other.

    A step has no downstream edge and so nothing to delimit; it is left out.
    """
    widths = sorted({w for _, _, _, w in entries if 0.0 < w < np.inf})
    if not widths:
        return
    if len(widths) == 1:
        ax.axvspan(*x_to_re(np.array([0.0, widths[0]])), color=SHADE_GREY,
                   lw=0, zorder=0)
        return
    col_of = {w: c for _, _, c, w in entries if 0.0 < w < np.inf}
    lo = 0.0
    for w in widths:
        ax.axvspan(*x_to_re(np.array([lo, w])), color=col_of[w],
                   alpha=SHADE_ALPHA, lw=0, zorder=0)
        lo = w


def x_to_re(x):
    return RE * np.sqrt(np.maximum(1.0 + BLASIUS_C**2 * np.asarray(x) / RE, 0.0))


def re_to_x(r):
    return ((np.asarray(r) / RE) ** 2 - 1.0) * RE / BLASIUS_C**2


def load_modes(name):
    path = os.path.join(DATA_ROOT, f"modes_Re{RE}_{name}_{TAG}.npz")
    if not os.path.isfile(path):
        return None
    z = np.load(path, allow_pickle=True)
    return ModeSet(
        x=z["x"],
        F=z["F"],
        omega=z["omega"],
        alpha=z["alpha"],
        n_global=int(z["n_global"]),
        meta=dict(z["meta"].item()),
    )


def fill_lost(F, g):
    """
    Fill the (station, F) pairs where the mode tracker returned nothing.

    ``track_modes`` writes NaN where it could not continue a branch. Across
    these fourteen cases that happens only on the F = 10 row of the three
    shallowest gaps -- one point at d = 0.25, twelve at d = 0.5, eighteen at
    d = 0.75 -- at the stations just downstream of the upstream edge where the
    lower branch is dropping out of the bottom of the scanned band and the
    mode becomes hard to follow.

    A NaN left in place cuts the column in two. The contour cannot close the
    lobe through it, the F = 0 row is written as though the column were stable
    there, and ``extrapolate_low_branch`` propagates it further still, since
    the value it places below the branch is computed from the row above it.
    The hole is filled by interpolating along F within the column, which is
    the step ``plots/localNfactorMaps.clean_map`` takes before any other.
    """
    g = np.array(g, dtype=float)
    for i in range(g.shape[1]):
        col = g[:, i]
        ok = np.isfinite(col)
        if ok.any() and not ok.all():
            col[~ok] = np.interp(F[~ok], F[ok], col[ok])
    return g


def despike(F, g, c):
    """
    -alpha_i with the holes of the mode selection filled, for contouring only.

    ``pp.localNfactor.select_mode`` rejects any eigenvalue damped by more than
    RATIO_MAX = 0.9 per radian of streamwise phase, and where the true mode is
    rejected it does not return nothing -- it returns the least damped *other*
    admissible eigenvalue, in practice the continuous spectrum at c_r ~ 0.73.
    That leaves weakly damped points sitting between two strongly amplified
    neighbours, which cut hairline slivers across the rise of the neutral
    curve at the upstream edge of the gap.

    Runs of such fallback points (recognised by c_r >= C_FALLBACK, not by
    their growth rate) that are bracketed in F by genuinely amplified modes
    are interpolated across, since what the scan reports there is that it
    found no instability mode, which is not the same as finding a stable one.
    Fallback points that are not so bracketed are left alone: away from the
    unstable region they are the honest answer.

    This is the ``despike`` step of ``plots/localNfactorMaps.clean_map`` --
    none of the neutral curve shaping the panels of that figure apply is done
    here.
    """
    g = np.array(g, dtype=float)
    fallback = ~(np.isfinite(c) & (c < C_FALLBACK))

    for i in range(g.shape[1]):
        j = 0
        while j < F.size:
            if not fallback[j, i]:
                j += 1
                continue
            k = j
            while k < F.size and fallback[k, i]:
                k += 1
            if j > 0 and k < F.size and g[j - 1, i] > 0 and g[k, i] > 0:
                g[j:k, i] = np.interp(F[j:k], [F[j - 1], F[k]], [g[j - 1, i], g[k, i]])
            j = k
    return g


def smooth_x(x, g, w, passes=None):
    """
    [1 2 1]/4 along x, applied within each streamwise segment separately.

    The segments are split at the gap edges, so the smoothing never reaches
    across x = 0 or x = w. Those two are real discontinuities of the base flow
    and the step the contour takes there is the answer, not an artifact; what
    this is for is the staircase left where the contour hands over from the
    shear layer branch to the Tollmien-Schlichting one, which happens well
    inside a segment.
    """
    passes = SMOOTH_PASSES if passes is None else passes
    if passes <= 0:
        return g
    out = np.array(g, dtype=float)
    segments = [x <= 0.0]
    if np.isfinite(w) and w > 0.0:
        segments += [(x > 0.0) & (x < w), x >= w]
    else:
        segments += [x > 0.0]
    for m in segments:
        idx = np.flatnonzero(m)
        if idx.size < 3:
            continue
        sm = out[:, idx]
        for _ in range(passes):
            p = np.pad(sm, ((0, 0), (1, 1)), mode="edge")
            sm = 0.25 * p[:, :-2] + 0.5 * p[:, 1:-1] + 0.25 * p[:, 2:]
        out[:, idx] = sm
    return out


def _low_branch(F, g, ing):
    """
    F of the lowest neutral crossing at each in gap station, NaN where the
    column is amplified all the way down to the bottom of the band and the
    branch therefore leaves through it rather than crossing.
    """
    out = np.full(g.shape[1], np.nan)
    for i in ing:
        col = g[:, i]
        j = np.flatnonzero(col[:-1] * col[1:] < 0.0)
        if j.size >= 2:
            k = j[0]
            out[i] = F[k] + (F[k + 1] - F[k]) * col[k] / (col[k] - col[k + 1])
    return out


def extrapolate_low_branch(x, F, g, w):
    """
    Continue the lower branch into the stations near the upstream edge where
    it drops below the scanned band and the scan loses it.

    A parabola is fitted to the first LOW_BRANCH_FIT stations that do carry
    the branch and matched in value at the junction, then the rows under the
    extrapolated branch are damped -- the row just below it by the amount that
    puts the crossing of the contour on the branch itself. The curve is then
    drawn by the contour like any other piece of the neutral curve, rather
    than laid on top of a field that still reads as amplified.

    This is a reading of the trend, not a computed neutral point.
    """
    ing = np.flatnonzero((x > 0.0) & (x < w))
    if not ing.size:
        return g
    lb = _low_branch(F, g, ing)
    have = ing[np.isfinite(lb[ing])]
    if have.size < 3:
        return g
    gone = ing[~np.isfinite(lb[ing])]
    gone = gone[gone < have.min()]
    if not gone.size:
        return g

    fit = have[:LOW_BRANCH_FIT]
    coef = np.polyfit(x[fit], lb[fit], 2)
    shift = lb[have[0]] - np.polyval(coef, x[have[0]])
    for i in gone:
        f_ex = np.polyval(coef, x[i]) + shift
        below = np.flatnonzero(F < f_ex)
        if not below.size or below.max() + 1 >= F.size:
            # the branch leaves through the bottom of the frame: leave the
            # column amplified all the way down and let it close on F = 0
            continue
        k = below.max()
        t = (f_ex - F[k]) / (F[k + 1] - F[k])
        g[k, i] = -t * g[k + 1, i] / (1.0 - t)
        g[:k, i] = -np.abs(g[:k, i])
    return g


def prepare(ms, name, w):
    """
    The field to contour, the F axis it sits on, and which stations close the
    lobe on F = 0.

    Order matters: the smoothing acts on the physical band alone, the F = 0
    row is added after it so the marker epsilon is not smeared into the data,
    and the extrapolation comes last because it needs that row to be able to
    place a branch below F = 10.
    """
    F = np.asarray(ms.F, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = ms.omega[:, None] / ms.alpha.real
    g = fill_lost(F, ms.growth)
    g = despike(F, g, c)
    g = smooth_x(ms.x, g, w)

    join = None
    if CLOSE_AT_F0:
        # the F = 0 row is neutral by construction; a signed epsilon carries
        # which side of the branch each station is on and puts the crossing of
        # the lower branch on the axis itself
        eps = 1e-9 * np.nanmax(np.abs(g))
        F = np.concatenate(([0.0], F))
        g = np.vstack((np.where(g[0] > 0.0, eps, -eps), g))

        if name in EXTRAPOLATE:
            g = extrapolate_low_branch(ms.x, F, g, w)
        join = g[0] > 0.0

    return F, g, join


def plot_family(key, fam, save_path):
    frame = FigureFrame(frame_w=4.7, aspect_ratio=0.75, pad_l=1.0, pad_b=2.1, pad_t=0.7)
    fig, ax = frame.fig, frame.ax
    ax.grid(False)

    ramp = plt.get_cmap(fam["cmap"])(np.linspace(0.38, 1.0, len(fam["cases"])))
    entries = [("flat", "Flat plate", COL_FLAT, gap_width("flat"))]
    entries += [(n, lbl, c, gap_width(n)) for n, lbl, c in fam["extras"]]
    entries += [
        (n, lbl, ramp[i], gap_width(n))
        for i, (n, lbl) in enumerate(fam["cases"].items())
    ]

    shade_gaps(ax, entries)

    missing = []
    for name, label, col, w in entries:
        ms = load_modes(name)
        if ms is None:
            missing.append(name)
            continue
        lw = LW_FLAT if name == "flat" else LW_CASE
        F, g, join = prepare(ms, name, w)
        r = x_to_re(ms.x)
        ax.contour(
            r,
            F,
            np.ma.masked_invalid(g),
            levels=[0.0],
            colors=[col],
            linewidths=lw,
            zorder=3 if name == "flat" else 2,
        )
        # the piece of the lobe that closes on F = 0 is the edge of the data,
        # not a zero crossing the contour can find, so it is drawn explicitly
        if join is not None and join.any():
            ax.plot(r, np.where(join, 0.0, np.nan), color=col, lw=lw,
                    zorder=3 if name == "flat" else 2)
        # contour draws no legend handle of its own
        ax.plot([], [], color=col, lw=lw, label=label)

    ax.set_xticks(RE_TICKS)
    ax.set_yticks(FTICKS)
    ax.set_xlim(*x_to_re(np.array(XLIM_X)))
    ax.set_ylim(*FLIM)
    ax.set_xlabel(r"$\mbox{\textit{Re}}_\ell$", labelpad=2)
    ax.set_ylabel(r"$F$", rotation=0, labelpad=8)

    # each horizontal edge carries its own ladder, major ticks only: Re_l on
    # the bottom, x on the top
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="both", top=False, right=True)
    ax.tick_params(axis="y", which="minor", length=2)
    frame.axis_on_top()

    ax2 = ax.secondary_xaxis("top", functions=(re_to_x, x_to_re))
    ax2.set_xlabel(r"$x$", labelpad=2.5)
    ax2.set_xticks(X_TICKS)
    ax2.xaxis.set_minor_locator(NullLocator())

    # the step is drawn early, so the gap curves land on top of it, but it is
    # listed last, after the family it is the limit of
    handles, labels = ax.get_legend_handles_labels()
    order = sorted(range(len(labels)), key=lambda i: labels[i] == "BFS")
    ax.legend(
        [handles[i] for i in order],
        [labels[i] for i in order],
        loc="upper center",
        bbox_to_anchor=(fam.get("legend_x", 0.5), -0.18),
        ncol=fam["legend_ncol"],
        columnspacing=1.0,
        frameon=False,
    )

    fig.savefig(save_path, format="pdf")
    plt.close(fig)
    print(f"Saved {os.path.normpath(save_path)}")
    if missing:
        print("  no scan for: " + ", ".join(missing))


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    out_dirs = [
        os.path.join(script_path, "../../../images"),
        os.path.join(script_path, "../../../latex/papers/jfm_incNS/Images"),
    ]
    import shutil

    for key, fam in FAMILIES.items():
        first = os.path.join(out_dirs[0], fam["out"])
        plot_family(key, fam, first)
        for d in out_dirs[1:]:
            shutil.copyfile(first, os.path.join(d, fam["out"]))


if __name__ == "__main__":
    main()
