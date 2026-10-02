"""
Local stability map of a gap flow: neutral curve, the growth rate it encloses
and the integrated n factor, one panel per (Reynolds number, geometry).

The growth rate is filled on a ladder of levels shared by every panel -- see
GROWTH_LEVELS -- so the colour alone says which Reynolds number is the more
unstable inside its own neutral curve. Its key is written once, horizontally,
to ``localMapGrowthColorbar.pdf``, to be placed wherever the manuscript wants
it.

Everything comes from the cached spatial Orr-Sommerfeld scans written by
``scripts/predictNfactorLocal.py`` (``data/localNfactor/modes_Re*_*.npz``),
which hold alpha(x, F) for a band of real reduced frequencies F = omega*1e6/Re.

Each panel carries four axes:

    bottom  Re      local displacement thickness Reynolds number
    top     x       streamwise station, measured from the upstream gap edge
    left    F       reduced frequency of the contour map
    right   n       n factor of the synthesised perturbation

so the (Re, F) stability diagram and the n(x) curve it integrates to share one
frame and can be read against each other station by station.
"""

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import AutoMinorLocator, FixedFormatter, NullLocator

from pp.figureFrame import FigureFrame
from pp.localNfactor import ModeSet, synthesise

plt.style.use("plots/style/jfm.mplstyle")
# the panels are laid out three across a 32pc JFM text block, so each PDF has
# to come out at exactly 4.6 cm *including* its labels: constrained layout at a
# fixed figure size, not a tight bounding box around an oversized canvas
plt.rcParams["savefig.bbox"] = "standard"

BLASIUS_C = 1.7207876573

REYNOLDS = [800, 1000, 3000]

DATA_ROOT = "/home/victor/Desktop/PhD/data/localNfactor"
TAG = "F2-600x40_n81_y45_X400"

# ---------------------------------------------------------------------------
# panel settings: one block per geometry, everything tweakable in one place
#
# The two geometries contrasted in the paper are a narrow deep gap and a wide
# shallow one, both on the stable side of the Hopf boundary at every Reynolds
# number considered. They need different frames: the deep gap carries a shear
# layer mode amplified up to the top of the scanned F band, the shallow one is
# neutral above F ~ 375, and the streamwise window of each is set by its own
# width (4w downstream of the upstream edge).
#
#   xlim / xticks   streamwise window and top axis ticks, in x
#   re_ticks        bottom axis ticks, in Re_l, per Reynolds number {re: [...]}.
#                   Each ladder starts at the Reynolds number of the run, so
#                   x = 0 on the top edge and Re_l = Re on the bottom one mark
#                   the same station, the upstream gap edge
#   flim / fticks   contour map, left axis. Both this and nlim start at zero:
#                   the two y axes share the bottom spine, so F = 0 and n = 0
#                   are the same line
#   nlim / nticks   n factor curves, right axis
#   nlim_re         per Reynolds number override of nlim, {re: (lo, hi)}, for
#   nticks_re       cases where one frame cannot hold every curve; nticks_re
#                   likewise. Both may be left empty
#
# The two horizontal ladders are kept apart: Re_l ticks on the bottom edge
# only, x ticks on the top edge only, major ticks alone on both.
#
# The last three switch the neutral curve clean up of ``clean_map``; read its
# docstring before turning any of them on for a new case, and in particular
# before trusting ``monotonic_gap``, which imposes a shape rather than
# reporting one:
#
#   close_at_f0     close the cavity branch on the F = 0 axis and mark it
#   f0_pad          how far below F = 0 the view reaches; 0 puts F = 0 on the
#                   bottom spine, so the branch marked on the axis is drawn on
#                   top of the spine rather than clear of it
#   despike         drop isolated (x, F) points where mode selection failed
#   monotonic_gap   force the in gap growth rate to be non decreasing in x
#   low_branch      per Reynolds number treatment of the lower branch inside
#                   the gap, "data" (the default) / "axis" / "extrapolate";
#                   see clean_map
#   low_branch_fit  stations used by the "extrapolate" fit
# ---------------------------------------------------------------------------
GEOMETRIES = {
    "d4_w12": dict(
        d=4.0,
        w=12.0,
        xlim=(-10.0, 4 * 12.0),
        xticks=[0, 12, 30, 45],
        re_ticks={
            800: [800, 825, 850],
            1000: [1000, 1025, 1050],
            3000: [3000, 3025, 3050],
        },
        flim=(0.0, 600.0),
        fticks=[0, 200, 400, 600],
        nlim=(0.0, 2.0),
        nticks=[0, 1, 2],
        nlim_re={},
        nticks_re={},
        close_at_f0=True,
        f0_pad=0.0,
        despike=True,
        monotonic_gap=True,
        low_branch={},
        low_branch_fit=20,
    ),
    "d0.5_w90": dict(
        d=0.5,
        w=90.0,
        xlim=(-70.0, 4 * 90.0),
        xticks=[0, 90, 200, 300],
        re_ticks={
            800: [800, 1000, 1200],
            1000: [1000, 1200, 1400],
            3000: [3000, 3200, 3400],
        },
        flim=(0.0, 450.0),
        fticks=[0, 200, 400],
        nlim=(0.0, 7.5),
        nticks=[0, 2, 4, 6],
        nlim_re={},
        nticks_re={},
        close_at_f0=True,
        f0_pad=0.0,
        despike=True,
        monotonic_gap=False,
        low_branch={800: "extrapolate", 1000: "extrapolate", 3000: "axis"},
        low_branch_fit=20,
    ),
}

# Filled levels of -alpha_i inside the neutral curve. The Tollmien-Schlichting
# growth rate outside the gap is O(5e-3) while the cavity shear layer reaches
# O(5e-1), so only a geometric ladder resolves both on one map: three decades,
# four bands each, boundaries landing on the decades themselves.
#
# The ladder is fixed rather than per panel, which is the whole point of the
# fill -- the same colour means the same growth rate in every panel, so the
# Re = 3000 lobe can be read against the Re = 800 one -- and it is *not*
# extended at the bottom, so everything below 1e-3 (in particular the whole
# stable side) is simply left unfilled.
GROWTH_LEVELS = np.geomspace(1e-3, 1.0, 13)
GROWTH_TICKS = [1e-3, 1e-2, 1e-1, 1.0]
GROWTH_TICKLABELS = [r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$", r"$10^{0}$"]
# light to dark, sequential: weak growth stays pale so the amplified core is
# what carries the ink, and no part of it competes with the red n factor curve
GROWTH_CMAP = "YlGnBu"
GROWTH_NORM = BoundaryNorm(GROWTH_LEVELS, plt.get_cmap(GROWTH_CMAP).N)

# the panels are laid out three across a 32pc JFM text block, so the PDF has to
# come out at exactly 4.5 x 3.9 cm *including* its labels: the frame and the
# paddings around it (all in cm) add up to that, and the frame is the same in
# every panel whatever its tick labels
FRAME_W_CM = 2.78
FRAME_ASPECT = 2.4 / FRAME_W_CM  # frame 2.4 cm high
PAD_L_CM = 0.86
PAD_R_CM = PAD_L_CM
PAD_B_CM = 0.75
PAD_T_CM = 0.65

# the standalone colour bar, to be placed wherever the manuscript wants it
CBAR_W_CM = 4.6
CBAR_H_CM = 1.1

COL_N = "tab:red"
COL_NEUTRAL = "black"
# the neutral curve is kept thin so that it reads as the edge of the filled
# region rather than as a separate object drawn over it
LW_NEUTRAL = 0.5


# phase speed at which the discretised continuous spectrum piles up. It is the
# C_MAX of pp.localNfactor.select_mode: when no genuine discrete mode is
# admissible the selection falls back on that pile up, so c_r >= this value
# means "no instability mode was found here", not "the mode runs this fast".
C_FALLBACK = 0.70


def clean_map(ms, geom, re):
    """
    Neutral curve clean up, applied to the contoured field only.

    Two of the three steps repair a known failure of the mode selection in
    ``pp.localNfactor.select_mode``, which rejects any eigenvalue damped by
    more than ``RATIO_MAX = 0.9`` per radian of streamwise phase. The cavity
    shear layer mode crosses that ratio as F -> 0, because alpha_r vanishes
    with the frequency while the growth rate does not: at Re = 800, x = 4.5 it
    is already 0.87 at F = 17 and reaches 1.4 at F = 2, where a direct solve
    with the cut lifted still returns a clean mode (alpha = 0.0066 - 0.0092i,
    c_r = 0.24, amplified). Where the true mode is rejected the selection does
    not return nothing, it returns the least damped *other* admissible
    eigenvalue -- in practice the continuous spectrum at c_r = 0.73 -- so the
    raw map is left with

      * NaN at the lowest frequencies, which clips the lower branch of the
        cavity lobe at the bottom of the grid instead of letting it close;
      * holes inside the amplified region, weakly damped points between two
        strongly amplified neighbours, which cut hairline slivers across the
        rise of the neutral curve at the upstream edge of the gap.

    ``despike`` fills the holes: runs of fallback points (recognised by
    c_r >= C_FALLBACK, not by their growth rate) that are bracketed in F by
    genuinely amplified modes are interpolated across, since what the scan
    actually reports there is that it found no instability mode, which is not
    the same as finding a stable one. Fallback points that are *not* so
    bracketed are left alone: away from the unstable region they are the
    honest answer, and at Re = 3000 they cover most of the map.

    ``close_at_f0`` closes the lobe on the axis. The lower branch really does
    sit at F = 0 inside the cavity: the growth rate falls smoothly to zero
    with the frequency and stays positive all the way down, which is the
    inflectional shear layer being unstable at every frequency it carries.

    ``low_branch`` decides what happens to the lower branch inside the gap,
    where the F band is too coarse to carry it: the scan puts its first point
    at F = 2 and its second at F = 17, and the branch lives in between. With
    "data" the contour is left as it comes out. "axis" flattens it onto F = 0
    across the whole gap, for cases where what survives of it is a sliver a
    few F units high -- at Re = 3000 the wide gap branch sits at F = 3.6 to
    4.4, one per cent of the panel, which is the axis to within the grid.
    "extrapolate" continues the branch into the stations near the upstream
    edge where it drops below F = 17 and the scan loses it, by fitting a
    parabola to the first ``low_branch_fit`` stations that do carry it and
    matching value at the junction. That extrapolation is a reading of the
    trend, not a computed neutral point, although re-solving those stations
    with the phase speed floor C_MIN lowered from 0.15 to 0.03 -- which is
    what loses them, the mode there running at c_r = 0.15 to 0.2 -- puts the
    branch within 1 to 3 units of F of it.

    ``monotonic_gap`` is a different kind of step and is *not* an artifact
    repair. The upper branch genuinely wanders inside the cavity -- at
    Re = 3000 it falls by 10 per cent over the upstream half before rising --
    and this replaces the in gap growth rate by its running maximum in x, so
    that the curve can only grow. What the panel then shows inside the grey
    band is an envelope of the local result, not the local result itself.

    The n factor curves are untouched: they are synthesised from the raw alpha
    in ``plot_panel``, so Delta N stays the number the scan produced.
    """
    x, F = ms.x, np.array(ms.F, dtype=float)
    g = np.array(ms.growth, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = ms.omega[:, None] / ms.alpha.real

    # frequencies the tracker lost, filled along F so the contour is not cut
    for i in range(g.shape[1]):
        col = g[:, i]
        ok = np.isfinite(col)
        if ok.any() and not ok.all():
            col[~ok] = np.interp(F[~ok], F[ok], col[ok])

    if geom["despike"]:
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
                bracketed = j > 0 and k < F.size and g[j - 1, i] > 0 and g[k, i] > 0
                if bracketed:
                    g[j:k, i] = np.interp(
                        F[j:k], [F[j - 1], F[k]], [g[j - 1, i], g[k, i]]
                    )
                j = k

    if geom["monotonic_gap"]:
        idx = np.flatnonzero((x > 0.0) & (x < geom["w"]))
        if idx.size:
            g[:, idx] = np.maximum.accumulate(g[:, idx], axis=1)

    eps = 1e-9 * np.nanmax(np.abs(g))
    mode = geom["low_branch"].get(re, "data")
    ing = np.flatnonzero((x > 0.0) & (x < geom["w"]))

    if mode == "axis" and ing.size:
        # everything under the upper branch is amplified, so the lobe is
        # bounded below by the axis and by nothing else
        for i in ing:
            pos = np.flatnonzero(g[:, i] > 0.0)
            if pos.size:
                g[: pos.max(), i] = np.maximum(g[: pos.max(), i], eps)

    elif mode == "extrapolate" and ing.size:
        lb = _low_branch(F, g, ing)
        have = ing[np.isfinite(lb[ing])]
        gone = ing[~np.isfinite(lb[ing])]
        gone = gone[gone < have.min()] if have.size else gone
        if have.size >= 3 and gone.size:
            fit = have[: geom["low_branch_fit"]]
            coef = np.polyfit(x[fit], lb[fit], 2)
            shift = lb[have[0]] - np.polyval(coef, x[have[0]])
            for i in gone:
                f_ex = np.polyval(coef, x[i]) + shift
                below = np.flatnonzero(F < f_ex)
                if not below.size or below.max() + 1 >= F.size:
                    continue  # the branch leaves through the axis: leave the
                    # column amplified all the way down
                k = below.max()
                # damp the rows under the extrapolated branch, the row just
                # below it by the amount that puts the crossing of the contour
                # on the branch itself. The curve is then drawn by the contour
                # like any other piece of the neutral curve, instead of being
                # laid on top of a field that still reads as amplified -- which
                # is what left a spike dropping to the axis at the junction.
                t = (f_ex - F[k]) / (F[k + 1] - F[k])
                g[k, i] = -t * g[k + 1, i] / (1.0 - t)
                g[:k, i] = -np.abs(g[:k, i])

    join = None
    if geom["close_at_f0"]:
        # the F = 0 row is neutral by construction; a signed epsilon carries
        # which side of the branch each station is on and puts the crossing of
        # the lower branch on the axis itself
        row0 = np.where(g[0] > 0.0, eps, -eps)
        F = np.concatenate(([0.0], F))
        g = np.vstack((row0, g))
        join = row0 > 0.0

    return F, g, join


def _low_branch(F, g, ing):
    """
    F of the lowest neutral crossing at each in gap station, NaN where the
    column is amplified all the way down to the bottom of the band and the
    branch therefore leaves through the axis rather than crossing.
    """
    out = np.full(g.shape[1], np.nan)
    for i in ing:
        col = g[:, i]
        j = np.flatnonzero(col[:-1] * col[1:] < 0.0)
        if j.size >= 2:
            k = j[0]
            out[i] = F[k] + (F[k + 1] - F[k]) * col[k] / (col[k] - col[k + 1])
    return out


def load_modes(re, name):
    path = os.path.join(DATA_ROOT, f"modes_Re{re}_{name}_{TAG}.npz")
    z = np.load(path, allow_pickle=True)
    # write data to file (not compressed) for inspection if needed
    

    return ModeSet(
        x=z["x"],
        F=z["F"],
        omega=z["omega"],
        alpha=z["alpha"],
        n_global=int(z["n_global"]),
        meta=dict(z["meta"].item()),
    )


def plot_panel(re, name, geom, save_path):
    ms = load_modes(re, name)
    ms_flat = load_modes(re, "flat")
    w = geom["w"]

    def x_to_re(x):
        return re * np.sqrt(np.maximum(1.0 + BLASIUS_C**2 * np.asarray(x) / re, 0.0))

    def re_to_x(r):
        return ((np.asarray(r) / re) ** 2 - 1.0) * re / BLASIUS_C**2

    frame = FigureFrame(
        frame_w=FRAME_W_CM,
        aspect_ratio=FRAME_ASPECT,
        pad_l=PAD_L_CM,
        pad_r=PAD_R_CM,
        pad_b=PAD_B_CM,
        pad_t=PAD_T_CM,
    )
    fig, ax = frame.fig, frame.ax
    ax.grid(False)

    r_grid = x_to_re(ms.x)
    F_map, g_map, join = clean_map(ms, geom, re)
    g = np.ma.masked_invalid(g_map)  # (nF, nx), -alpha_i

    # the gap itself, where the parallel flow assumption is at its weakest and
    # where essentially all of the excess amplification is produced
    ax.axvspan(x_to_re(0.0), x_to_re(w), color="0.85", lw=0, zorder=0)

    # growth rate inside the neutral curve. Not extended below the lowest
    # level, so the stable side stays unfilled and the edge of the fill is the
    # 1e-3 isoline, which hugs the neutral curve everywhere the growth rate
    # rises as fast as it does here
    cf = ax.contourf(
        r_grid,
        F_map,
        g,
        levels=GROWTH_LEVELS,
        norm=GROWTH_NORM,
        cmap=GROWTH_CMAP,
        extend="neither",
        zorder=1,
    )
    # without this the PDF shows white seams between bands, where the viewer
    # antialiases each filled polygon against its neighbour
    cf.set_edgecolor("face")

    # neutral curve
    ax.contour(
        r_grid,
        F_map,
        g,
        levels=[0.0],
        colors=COL_NEUTRAL,
        linewidths=LW_NEUTRAL,
        zorder=3,
    )
    # the branch that closes on the axis: drawn explicitly, since it is the
    # edge of the data and not a zero crossing the contour can find
    if join is not None and join.any():
        ax.plot(
            r_grid,
            np.where(join, 0.0, np.nan),
            color=COL_NEUTRAL,
            lw=LW_NEUTRAL,
            zorder=3,
        )

    xlim = geom["xlim"]
    ax.set_xlim(*x_to_re(np.array(xlim)))
    ax.set_xlabel(r"$\mbox{\textit{Re}}_\ell$", labelpad=2)
    ax.set_ylabel(r"$F$", rotation=0, labelpad=5)
    # set_yticks widens the view to cover the ticks, so the limits have to be
    # imposed after them, not before -- otherwise a tick above the requested
    # top silently stretches the axis (see the same order on axN below)
    ax.set_yticks(geom["fticks"])
    flo, fhi = geom["flim"]
    ax.set_ylim(flo - (geom["f0_pad"] if geom["close_at_f0"] else 0.0), fhi)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    # each horizontal edge carries its own ladder, major ticks only: Re_l on
    # the bottom, x on the top. The right edge belongs to the n axis, so ax
    # leaves it alone
    ax.set_xticks(geom["re_ticks"][re])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.tick_params(which="both", top=False, right=False)
    ax.tick_params(axis="y", which="minor", length=2)
    frame.axis_on_top()

    ax2 = ax.secondary_xaxis("top", functions=(re_to_x, x_to_re))
    ax2.set_xlabel(r"$x$", labelpad=2.5)
    ax2.set_xticks(geom["xticks"])
    ax2.xaxis.set_minor_locator(NullLocator())

    # n factor, synthesised from the same alpha_i by summing the modes in
    # energy with flat weights (see pp.localNfactor)
    axN = ax.twinx()
    axN.grid(False)
    xf, nf = synthesise(ms_flat, None)
    xg, ng = synthesise(ms, None)
    axN.plot(x_to_re(xf), nf, color=COL_N, lw=0.8, ls="--", zorder=4)
    axN.plot(x_to_re(xg), ng, color=COL_N, lw=1.2, ls="-", zorder=5)
    axN.set_yticks(geom["nticks_re"].get(re, geom["nticks"]))
    axN.set_ylim(*geom["nlim_re"].get(re, geom["nlim"]))
    axN.set_ylabel(r"$n$", rotation=0, labelpad=5, color=COL_N)
    axN.yaxis.set_minor_locator(AutoMinorLocator(2))
    axN.tick_params(axis="y", which="both", colors=COL_N, direction="in")
    axN.tick_params(axis="y", which="minor", length=2)
    for s in ("right",):
        axN.spines[s].set_color(COL_N)

    fig.savefig(save_path, format="pdf")
    plt.close(fig)
    msg = f"Saved {os.path.normpath(save_path)}  (dN_end = {ng[-1] - nf[-1]:.2f})"
    # only the part of each curve that is actually inside the streamwise window
    # can run off the top, so the warning is windowed too
    nhi = geom["nlim_re"].get(re, geom["nlim"])[1]
    peak = max(
        n[(x >= xlim[0]) & (x <= xlim[1])].max(initial=-np.inf)
        for x, n in ((xf, nf), (xg, ng))
    )
    if peak > nhi:
        msg += f"  [n clipped: curves reach {peak:.2f} > {nhi:g} in view]"
    print(msg)


def plot_colorbar(save_path):
    """The growth rate colour bar on its own, horizontal, for free placement."""
    fig, ax = plt.subplots(
        figsize=(CBAR_W_CM / 2.54, CBAR_H_CM / 2.54),
        layout="constrained",
    )
    fig.get_layout_engine().set(w_pad=0.015, h_pad=0.01, wspace=0, hspace=0)
    # the style sheet turns the grid on for every axes, and on a colour bar
    # that draws grey rules straight across the bands
    ax.grid(False)

    cb = fig.colorbar(
        ScalarMappable(norm=GROWTH_NORM, cmap=GROWTH_CMAP),
        cax=ax,
        orientation="horizontal",
        # the bands are geometric, so equal widths on the bar; "proportional"
        # would crush the three lower decades against the left end
        spacing="uniform",
        ticks=GROWTH_TICKS,
    )
    cb.set_label(r"$-\alpha_i$", labelpad=2)
    cb.ax.xaxis.set_major_formatter(FixedFormatter(GROWTH_TICKLABELS))
    cb.ax.tick_params(which="both", direction="out", length=2.0, width=0.5, pad=1.5)
    cb.outline.set_linewidth(0.5)
    # the map draws no band edges, so neither does the bar
    cb.dividers.set_visible(False)

    fig.savefig(save_path, format="pdf")
    plt.close(fig)
    print(f"Saved {os.path.normpath(save_path)}")


def main():
    import shutil

    script_path = os.path.dirname(os.path.abspath(__file__))
    out_dirs = [
        os.path.join(script_path, "../../../images"),
        os.path.join(script_path, "../../../latex/papers/jfm_incNS/Images"),
    ]

    def publish(fname, draw):
        first = os.path.join(out_dirs[0], fname)
        draw(first)
        for d in out_dirs[1:]:
            shutil.copyfile(first, os.path.join(d, fname))

    for re in REYNOLDS:
        for name, geom in GEOMETRIES.items():
            publish(
                f"localMapRe{re}_{name}.pdf",
                lambda p, re=re, name=name, geom=geom: plot_panel(re, name, geom, p),
            )
    publish("localMapGrowthColorbar.pdf", plot_colorbar)


if __name__ == "__main__":
    main()
