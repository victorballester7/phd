"""
N factor curves of the gap flow: the global result against the local
(parallel flow) prediction of the same case, one panel per family.

This is the companion of ``plots/nfactorCurves.py``. That figure carries the
global curves alone -- n(x) read off the linearised Navier-Stokes solution
forced by blowing and suction (``pp.DeltaN_computation.computeNx`` applied to
``data/pointsavg_n600.dat``). Here every case is drawn twice:

    solid   global      n from the direct linear solver
    dashed  local       n predicted by local stability theory

The local curve is the incoherent sum over the scanned frequency band,
``pp.localNfactor.synthesise`` with flat weights, taken from the cached
Orr-Sommerfeld scans of ``scripts/predictNfactorLocal.py``. It is the
like for like counterpart of the global curve: both are 1/2 log(A^2/A0^2) of
an amplitude built from every frequency at once, referenced to the same
upstream station, rather than the envelope max_j N_j of a single mode.

Two families, both at Re = 1000:

    d1.5    constant depth, w = 10 ... 50, plus the backward facing step of
            the same depth -- the w -> infinity limit the family walks towards
    w40     constant width, d = 0.25 ... 1.75

The flat plate is the reference both families are measured against and is
carried in both panels.

The bottom axis is the local displacement thickness Reynolds number
Re_l(x) = Re sqrt(1 + C^2 x / Re), the top axis the streamwise station x
measured from the upstream gap edge.
"""

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, NullLocator

from pp.DeltaN_computation import computeNx
from pp.fileManagement import extract_depth_width
from pp.localNfactor import ModeSet, synthesise
from pp.figureFrame import FigureFrame

plt.style.use("plots/style/jfm.mplstyle")

RE = 1000
BLASIUS_C = 1.7207876573

ROOT = "/home/victor/Desktop/PhD"
SRC_ROOT = os.path.join(ROOT, "src")
DATA_ROOT = os.path.join(ROOT, "data", "localNfactor")
NSAMPLE = 600

# the cached scan these curves are read from. It has to name a band, because
# the local n factor is a sum over the scanned frequencies and so depends on
# which ones were scanned: F = 10 to 500 in 28 steps, 81 collocation nodes,
# free stream at y = 45, stations out to x = 700.
TAG = "F10-500x28_n81_y45_X700"

# ---------------------------------------------------------------------------
# the two families
#
#   cases       cache / case name -> legend label
#   cmap        colour ramp walked across the family. It starts at 0.38, not
#               at 0, because below that the first case is too pale to survive
#               printing (the same floor plots/nfactorCurves.py uses)
#   extras      cases drawn in their own colour rather than off the ramp
#   out         file name of the panel
# ---------------------------------------------------------------------------
COL_FLAT = "tab:blue"
COL_BFS = "tab:orange"

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
        out="nfactorCurvesLocald1.5.pdf",
        legend_ncol=3,
    ),
    "w40": dict(
        cases={
            "d0.25_w40": r"$d=0.25$",
            "d0.5_w40": r"$d=0.50$",
            "d0.75_w40": r"$d=0.75$",
            "d1_w40": r"$d=1.00$",
            "d1.25_w40": r"$d=1.25$",
            "d1.5_w40": r"$d=1.50$",
            # "d1.75_w40_old": r"$d=1.75$",
            "d1.75_w40": r"$d=1.75$",
        },
        cmap="Reds",
        extras=[],
        out="nfactorCurvesLocalw40.pdf",
        legend_ncol=3,
    ),
}

# the same frame in both panels, so that (a) and (b) can be read against each
# other -- and the same one plots/nfactorCurves.py uses, so that this figure
# can be laid over its global only counterpart
# the top is 7.8 rather than the 7.2 of nfactorCurves.py: the local curve of
# the step reaches 7.31 at the downstream end, and the tick ladder is kept the
# same so the two figures still read as one pair
YLIM = (0.0, 7.5)
YTICKS = [0, 2, 4, 6]
XLIM_X = (-80.0, 600.0)  # in x; the Re_l window follows from it
# The two horizontal ladders are kept apart: Re_l ticks on the bottom edge
# only, x ticks on the top edge only, major ticks alone on both. The Re_l
# ladder starts at Re and the x one at 0, so the first tick of each marks the
# same station, the upstream gap edge.
RE_TICKS = [1000, 1200, 1400, 1600]
X_TICKS = [-70, 0, 40, 200, 400, 600]

LW_GLOBAL = 1.0
LW_LOCAL = 0.9
LS_GLOBAL = "-"
LS_LOCAL = "--"

# The window the n factor is averaged over, ticked on the global curve of each
# gap case as plots/nfactorCurves.py did: AVG_OFFSET to AVG_OFFSET + AVG_LENGTH
# past the downstream gap edge, in a darkened shade of the case colour. The
# ticks are AVG_NMARK stations spaced evenly in x across the window, with n
# interpolated onto them, rather than the raw data points, which are unevenly
# spaced and bunch up on some of the curves.
AVG_OFFSET = 80.0
AVG_LENGTH = 50.0
AVG_NMARK = 6
AVG_DARKEN = 0.85

# gap shading, drawn behind the curves
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

    One grey band when every case in the family has the same gap. When the
    widths differ, each case colours the stretch its own gap adds to the one
    before it -- the widths are nested, all starting at x = 0 -- so the
    shading is a staircase whose steps are the downstream edges, and each
    colour appears once, at its own edge, with no translucent bands stacked on
    top of each other.

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
    """Local displacement thickness Reynolds number of the Blasius station."""
    return RE * np.sqrt(np.maximum(1.0 + BLASIUS_C**2 * np.asarray(x) / RE, 0.0))


def re_to_x(r):
    return ((np.asarray(r) / RE) ** 2 - 1.0) * RE / BLASIUS_C**2


# ---------------------------------------------------------------------------
# the two n factors of a case
# ---------------------------------------------------------------------------
def _avg_file(name):
    """Averaged field of the global run, or None when it was never extracted."""
    if name == "flat":
        p = os.path.join(SRC_ROOT, f"flatPlateRe{RE}inc", "directLinearSolver",
                         "blowingSuction", "data", f"pointsavg_n{NSAMPLE}.dat")
    elif name.startswith("bfs_"):
        p = os.path.join(SRC_ROOT, f"bfsRe{RE}inc", "directLinearSolver",
                         "blowingSuction", name[4:], "data",
                         f"pointsavg_n{NSAMPLE}.dat")
    else:
        p = os.path.join(SRC_ROOT, f"incGapRe{RE}", "directLinearSolver",
                         "blowingSuction", name, "data",
                         f"pointsavg_n{NSAMPLE}.dat")
    return p if os.path.isfile(p) else None


def global_curve(name):
    """n(x) from the direct linear solver, or None when the run has no field."""
    path = _avg_file(name)
    if path is None:
        return None
    x, n = computeNx(path, False)
    if not len(x):
        return None
    if name == "flat":
        # the last station of the flat plate run sits where the 2D
        # extrapolation of the averaged field breaks down; it is replaced by
        # the linear continuation of the two before it, as in nfactorCurves.py
        n[-1] = 2.0 * n[-2] - n[-3]
    return x, n


def local_curve(name):
    """n(x) predicted by local theory, or None when the case was not scanned."""
    path = os.path.join(DATA_ROOT, f"modes_Re{RE}_{name}_{TAG}.npz")
    if not os.path.isfile(path):
        return None
    z = np.load(path, allow_pickle=True)
    ms = ModeSet(
        x=z["x"],
        F=z["F"],
        omega=z["omega"],
        alpha=z["alpha"],
        n_global=int(z["n_global"]),
        meta=dict(z["meta"].item()),
    )
    return synthesise(ms, None)


def plot_family(key, fam, save_path):
    frame = FigureFrame(frame_w=5.3, aspect_ratio=0.75, pad_l=0.71, pad_b=2.7, pad_t=0.7)
    fig, ax = frame.fig, frame.ax
    ax.grid(False)

    ramp = plt.get_cmap(fam["cmap"])(np.linspace(0.3, 1.0, len(fam["cases"])))

    # flat plate first, so the gap curves are drawn over it
    entries = [("flat", "Flat plate", COL_FLAT, gap_width("flat"))]
    entries += [(n, lbl, c, gap_width(n)) for n, lbl, c in fam["extras"]]
    entries += [
        (n, lbl, ramp[i], gap_width(n))
        for i, (n, lbl) in enumerate(fam["cases"].items())
    ]

    shade_gaps(ax, entries)

    missing = []
    for name, label, col, w in entries:
        g = global_curve(name)
        loc = local_curve(name)
        if g is None and loc is None:
            missing.append(f"{name} (both)")
            continue
        if g is None:
            missing.append(f"{name} (global)")
        if loc is None:
            missing.append(f"{name} (local)")

        if g is not None:
            ax.plot(x_to_re(g[0]), g[1], ls=LS_GLOBAL, lw=LW_GLOBAL, color=col,
                    label=label)
            if name in fam["cases"]:
                x0 = w + AVG_OFFSET
                xm = np.linspace(x0, x0 + AVG_LENGTH, AVG_NMARK)
                ax.plot(x_to_re(xm), np.interp(xm, g[0], g[1]), ls="", marker="|",
                        markersize=3, color=np.asarray(col) * AVG_DARKEN)
        if loc is not None:
            # no label: the dashed line shares its colour with the solid one,
            # and the linestyle is keyed by the second legend instead
            ax.plot(x_to_re(loc[0]), loc[1], ls=LS_LOCAL, lw=LW_LOCAL, color=col,
                    label=None if g is not None else label)

    ax.set_xticks(RE_TICKS)
    ax.set_yticks(YTICKS)
    ax.set_xlim(*x_to_re(np.array(XLIM_X)))
    ax.set_ylim(*YLIM)
    ax.set_xlabel(r"$\mbox{\textit{Re}}_\ell$", labelpad=2)
    ax.set_ylabel(r"$n$", rotation=0, labelpad=10)

    # each horizontal edge carries its own ladder, major ticks only: Re_l on
    # the bottom, x on the top
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="both", top=False, right=True)
    ax.tick_params(axis="y", which="minor", length=2)
    frame.axis_on_top()

    ax2 = ax.secondary_xaxis("top", functions=(re_to_x, x_to_re))
    ax2.set_xlabel(r"$x$", labelpad=3)
    ax2.set_xticks(X_TICKS)
    ax2.xaxis.set_minor_locator(NullLocator())

    # Two legends, both under the axes: colour says which case, linestyle says
    # which analysis. Neither corner inside the frame is free -- the curves
    # climb to the top right and the flat plate runs along the bottom right --
    # so a legend placed inside would sit on data in one family or the other.
    style_legend = ax.legend(
        handles=[
            Line2D([], [], color="0.25", ls=LS_GLOBAL, lw=LW_GLOBAL, label="Global LSA"),
            Line2D([], [], color="0.25", ls=LS_LOCAL, lw=LW_LOCAL, label="Local LSA"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=2,
        # columnspacing=1.0,
        # handlelength=1.8,
        # handletextpad=0.5,
        frameon=False,
    )
    ax.add_artist(style_legend)
    # the step is drawn early, so the gap curves land on top of it, but it is
    # listed last, after the family it is the limit of
    handles, labels = ax.get_legend_handles_labels()
    order = sorted(range(len(labels)), key=lambda i: labels[i] == "BFS")
    ax.legend(
        [handles[i] for i in order],
        [labels[i] for i in order],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.3),
        ncol=fam["legend_ncol"],
        columnspacing=1.0,
        # handlelength=1.6,
        # handletextpad=0.5,
        frameon=False
    )

    fig.savefig(save_path, format="pdf")
    plt.close(fig)
    print(f"Saved {os.path.normpath(save_path)}")
    if missing:
        print("  no curve for: " + ", ".join(missing))


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
