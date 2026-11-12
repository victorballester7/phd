"""
Overview of the (w, d) plane at Re = 1000 for the introduction.

A stripped-down version of 3d2dStabilityCurveRe1000.py. The reader meets this
figure before any notation is introduced, so the regions are shown as flat
tints labelled in words, not as numbered hatched regions, and the
interpolation bands and the zoom panel are left out.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Polygon
from pp.colors import colors, mix_with_black
from pp.hopfDataPoints import (
    hopf2d_re_stable,
    hopf_re_unstable,
    instability3d_stable,
    instability3d_unstable,
    exp_bypass_curve,
)
from pp.filterData import average_curves
from pp.figureFrame import FigureFrame


plt.style.use("plots/style/jfm.mplstyle")

COL_2D = "tab:green"
COL_3D = "tab:orange"
COL_EXP = "k"  # same as the 2D regime map, where brown was hard to pick out
COL_BOTH = "tab:blue"  # neither tint: blending green and orange reads as mud
ALPHA_FILL = 0.2
N_PTS = 2000

W_LIM = (0, 130)
D_LIM = (0, 4)

# Extra space between the frame and the legend and gap sketch on its right
# (cm). Sized so that, once the figure is centred in the JFM text width, the
# legend just right of the frame; gap sketch in the left pad, its right edge
# this far from the frame (clear of the y tick labels)
LEGEND_GAP = 0.0  # cm
SKETCH_GAP = 0.45  # cm


def _curve_y(x, curve):
    """d(w) of an (N, 2) boundary. Left of the curve's first point no gap is
    unstable, so the boundary is sent above the frame."""
    order = np.argsort(curve[:, 0])
    return np.interp(x, curve[order, 0], curve[order, 1], left=2 * D_LIM[1])


def fillRegions(ax, inst3d, hopf2d):
    """Tint the three unstable regions: 3D only, 2D only, and both."""
    x = np.linspace(*W_LIM, N_PTS)
    y2d = np.minimum(_curve_y(x, hopf2d), D_LIM[1])
    y3d = np.minimum(_curve_y(x, inst3d), D_LIM[1])
    lo, hi = np.minimum(y2d, y3d), np.maximum(y2d, y3d)
    kw = dict(alpha=ALPHA_FILL, linewidth=0, zorder=1)
    ax.fill_between(x, lo, hi, where=y3d <= y2d, facecolor=COL_3D, interpolate=True, **kw)
    ax.fill_between(x, lo, hi, where=y2d <= y3d, facecolor=COL_2D, interpolate=True, **kw)
    ax.fill_between(x, hi, D_LIM[1], facecolor=COL_BOTH, **kw)


def plotCurves(ax, inst3d, hopf2d):
    ax.plot(hopf2d[:, 0], hopf2d[:, 1], "-", color=COL_2D, lw=1.3, zorder=3)
    ax.plot(inst3d[:, 0], inst3d[:, 1], "-.", color=COL_3D, lw=1.3, zorder=3)
    ax.plot(
        exp_bypass_curve[:, 0],
        exp_bypass_curve[:, 1],
        "--",
        color=COL_EXP,
        lw=1.3,
        zorder=4,
    )


def addText(ax):
    """Region labels in plain words; the two thin strips get an arrow."""
    k = 0.45
    leader = lambda c: dict(
        arrowstyle="-|>,head_length=0.4,head_width=0.18",
        color=mix_with_black(c, k),
        lw=0.6,
        shrinkA=2,
        shrinkB=0,
    )

    ax.text(
        70,
        0.45,
        "steady flow: the gap only\n\\textbf{amplifies} incoming \\textbf{TS waves}",
        ha="center",
        va="center",
        color="0.25",
    )
    ax.text(
        80,
        3.0,
        "unstable to 2D and 3D perturbations:\n\\textbf{transition at the gap}",
        ha="center",
        va="center",
        color=mix_with_black(COL_BOTH, k),
    )
    ax.annotate(
        "3D only:\n\\textbf{centrifugal}\n\\textbf{instability}",
        xy=(15, 2.6),
        xytext=(3, 1.0),
        color=mix_with_black(COL_3D, k),
        ha="left",
        va="center",
        arrowprops=leader(COL_3D),
    )
    ax.annotate(
        "2D only:\n\\textbf{Rossiter}-like\n\\textbf{oscillation}",
        xy=(94, 1.65),
        xytext=(127, 0.95),
        color=mix_with_black(COL_2D, k),
        ha="right",
        va="center",
        arrowprops=leader(COL_2D) | dict(shrinkA=1, relpos=(0.35, 1.0)),
    )


def addLegend(frame):
    from matplotlib.lines import Line2D

    handles = [
        Line2D([0], [0], color=COL_2D, ls="-", lw=1.3),
        Line2D([0], [0], color=COL_3D, ls="-.", lw=1.3),
        Line2D([0], [0], color=COL_EXP, ls="--", lw=1.3),
    ]
    labels = [
        "2D stability\nboundary",
        "3D stability\nboundary",
        "Experimental\nbypass\ntransition\nboundary\n(Crouch\n$\\it{et\\ al.}$ 2022)",
    ]
    frame.ax.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(0.985 + LEGEND_GAP / frame.frame_w, 0.5),
        frameon=False,
        handlelength=2.2,
        handletextpad=0.5,
    )


def addGapSketch(frame):
    """Small rectangular gap defining w, d and dstar, in the left pad at the
    top of the frame."""
    xl, xr, yb, yt = -1.5, 2.3, -1.08, 1.13  # data limits
    W, H = (xr - xl) / 1.5, (yt - yb) / 1.7  # cm
    sk = frame.add_axes_cm(
        frame.pad_l - SKETCH_GAP - W,
        frame.pad_b + frame.frame_h + 0.05 - H,  # top-aligned with the frame
        W,
        H,
    )
    sk.set_axis_off()
    sk.set_xlim(xl, xr)
    sk.set_ylim(yb, yt)
    gw, gd = 1.5, 0.9  # gap width and depth
    wall = [(xl, 0), (0, 0), (0, -gd), (gw, -gd), (gw, 0), (xr, 0)]
    sk.add_patch(
        Polygon(wall + [(xr, yb), (xl, yb)], closed=True, fc="0.85", ec="none")
    )
    sk.plot(*zip(*wall), color="k", lw=0.8)
    # displacement thickness growing up to the upstream edge, where it is marked
    xs = np.linspace(xl, 0, 50)
    ds = 0.5
    sk.plot(xs, ds * np.sqrt((xs + 3) / 3), color="k", lw=0.5, ls=(0, (2, 1.5)))
    sk.annotate(
        "", xy=(-0.55, 0), xytext=(-0.55, ds),
        arrowprops=dict(arrowstyle="<->,head_length=0.25,head_width=0.12", color="k", lw=0.5, shrinkA=0, shrinkB=0),
    )
    sk.text(-0.65, 0.2, r"$\delta^*$", ha="right", va="center")
    # boundary-layer flow arrow with a stealth head (matplotlib has no built-in
    # one); head sized in cm, converted to data units since x and y scale differ
    x0, x1, y = -1.4, -0.3, 1.05
    cx, cy = W / (xr - xl), H / (yt - yb)  # cm per data unit
    hl, hw = 0.2 / cx, 0.07 / cy  # head length, half-width
    sk.plot([x0, x1 - 0.5 * hl], [y, y], color="k", lw=0.6)
    sk.add_patch(
        Polygon(
            [(x1, y), (x1 - hl, y + hw), (x1 - 0.6 * hl, y), (x1 - hl, y - hw)],
            closed=True, fc="k", ec="none",
        )
    )
    # dimension marks
    sk.annotate(
        "", xy=(0, 0.15), xytext=(gw, 0.15),
        arrowprops=dict(arrowstyle="<->", color="k", lw=0.5, shrinkA=0, shrinkB=0),
    )
    sk.text(gw / 2, 0.22, r"$w$", ha="center", va="bottom")
    sk.annotate(
        "", xy=(gw + 0.2, 0), xytext=(gw + 0.2, -gd),
        arrowprops=dict(arrowstyle="<->", color="k", lw=0.5, shrinkA=0, shrinkB=0),
    )
    sk.text(gw + 0.4, -gd / 2, r"$d$", ha="left", va="center")


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(script_path, "../../../images/introRegimeMap.pdf")

    frame = FigureFrame(
        frame_w=7.4, aspect_ratio=0.633, pad_l=3.05, pad_b=0.8, pad_t=0.11  # pad_r = pad_l: frame centred, legend and sketch in pad_r
    )
    fig, ax = frame.fig, frame.ax

    hopf2d = average_curves(hopf2d_re_stable[1000], hopf_re_unstable[1000], npts=500)[2]
    inst3d = average_curves(
        instability3d_stable[1000], instability3d_unstable[1000], npts=10000
    )[2]

    ax.set_xticks([0, 20, 40, 60, 80, 100, 120])
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_xlim(W_LIM)
    ax.set_ylim(D_LIM)
    ax.set_xlabel(r"$w$")
    ax.set_ylabel(r"$d$", rotation=0, labelpad=8)
    ax.grid(False)
    frame.jfm_ticks()
    frame.axis_on_top()

    fillRegions(ax, inst3d, hopf2d)
    plotCurves(ax, inst3d, hopf2d)
    addText(ax)
    addLegend(frame)
    addGapSketch(frame)

    fig.savefig(save_path, format="pdf")
    print(colors.OKGREEN + f"✓ Plot saved to: {save_path}" + colors.ENDC)


if __name__ == "__main__":
    main()
