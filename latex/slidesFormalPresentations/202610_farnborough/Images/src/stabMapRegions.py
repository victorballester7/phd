"""Slide thumbnails of the 2D/3D stability map (Re = 1000), each with some of
the regions I-V highlighted: no zoom inset, no legend, chosen regions filled.

Regions, curves and labels come from scripts/pp/plots/3d2dStabilityCurveRe1000.py
(and its NoInset variant), so the map stays identical to the paper's. Run from
scripts/pp:

    uv run python ~/Desktop/PhD/latex/slidesFormalPresentations/202610_farnborough/Images/src/stabMapRegions.py

Axes position in the output, as fractions of the image (for TikZ markers):
    x = (0.9 + 5.5*w/130)/6.7,  y = (0.8 + 4.125*d/4)/5.035
"""

import importlib.util
import os

import numpy as np

from pp.hopfDataPoints import (
    hopf2d_re_stable,
    hopf_re_unstable,
    instability3d_stable,
    instability3d_unstable,
    exp_bypass_curve,
)
from pp.filterData import average_curves
from pp.figureFrame import FigureFrame

PLOTS = os.path.join(os.getcwd(), "plots")
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

HIGHLIGHT = "tab:red"
HIGHLIGHT_ALPHA = 0.35

# output file -> regions to highlight
MAPS = {
    "stabMapRegionII.pdf": ["II"],
    "stabMapRegionIII_IV.pdf": ["III", "IV"],
    "stabMapRegionV.pdf": ["V"],
}


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(PLOTS, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


base = load("stab3d2d", "3d2dStabilityCurveRe1000.py")
noinset = load("stab3d2dNoInset", "3d2dStabilityCurveRe1000NoInset.py")


def regions(ax, inst3ds, hopf2ds):
    """The five regions exactly as outlined in base.makePatterns."""
    inst3d, hopf2d = inst3ds[2], hopf2ds[2]
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xs = np.linspace(0.0, np.max(inst3d[:, 0]), 2000)
    gap = base._curve_y(xs, inst3d) - base._curve_y(xs, hopf2d)
    x_cross = xs[np.argmin(np.abs(gap))]

    def bounds(xa, xb):
        x = np.linspace(xa, xb, base.N_REGION_PTS)
        y2d, y3d = base._curve_y(x, hopf2d), base._curve_y(x, inst3d)
        return x, np.minimum(y2d, y3d), np.maximum(y2d, y3d)

    x, lo, hi = bounds(x0, x1)
    return {
        "I": (x, np.full_like(x, y0), lo),
        "II": bounds(x0, x_cross),
        "III": bounds(x_cross, base.X_SEP_III_IV),
        "IV": bounds(base.X_SEP_III_IV, x1),
        "V": (x, hi, np.full_like(x, y1)),
    }


def plot_map(highlight, save_path):
    frame = FigureFrame(
        frame_w=5.5, aspect_ratio=0.75, pad_l=0.9, pad_r=0.3, pad_b=0.8, pad_t=0.11
    )
    fig, ax = frame.fig, frame.ax

    hopf2ds = average_curves(hopf2d_re_stable[1000], hopf_re_unstable[1000], npts=500)
    inst3ds = average_curves(
        instability3d_stable[1000], instability3d_unstable[1000], npts=10000
    )

    ax.set_xticks([0, 20, 40, 60, 80, 100, 120])
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_xlim([0, 130])
    ax.set_ylim([0, 4])
    ax.set_xlabel(r"$w$")
    ax.set_ylabel(r"$d$", rotation=0, labelpad=10)
    ax.grid(False)
    frame.jfm_ticks()
    frame.axis_on_top()

    base.makePatterns(ax, inst3ds, hopf2ds)
    for name, (x, ya, yb) in regions(ax, inst3ds, hopf2ds).items():
        if name in highlight:
            ax.fill_between(
                x, ya, yb,
                facecolor=HIGHLIGHT, alpha=HIGHLIGHT_ALPHA, linewidth=0, zorder=1.5,
            )
    base.plotCurves(ax, inst3ds, hopf2ds, exp_bypass_curve)
    noinset.addText(ax)

    fig.savefig(save_path, format="pdf")
    print("Saved", os.path.normpath(save_path))


if __name__ == "__main__":
    for fname, highlight in MAPS.items():
        plot_map(highlight, os.path.join(OUT_DIR, fname))
