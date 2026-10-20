"""Two overlay steps of 3d2dStabilityCurveRe1000.pdf for the slides, with
identical geometry so the TikZ markers drawn on top stay put:

    stabCurveStep1.pdf   whole map, zoom panel/box/connectors hidden
    stabCurveStep2.pdf   whole map with the zoom, box and connectors coloured

main() of scripts/pp/plots/3d2dStabilityCurveRe1000.py, with the zoom made
optional and recoloured. Run from scripts/pp:

    uv run python ~/Desktop/PhD/latex/slidesFormalPresentations/202610_farnborough/Images/src/stabCurveInsetSteps.py
"""

import importlib.util
import os

from matplotlib.transforms import TransformedBbox
from mpl_toolkits.axes_grid1.inset_locator import BboxConnector, BboxPatch

from pp.hopfDataPoints import (
    hopf2d_re_stable,
    hopf_re_unstable,
    instability3d_stable,
    instability3d_unstable,
    exp_bypass_curve,
)
from pp.filterData import average_curves
from pp.figureFrame import FigureFrame, axis_on_top, jfm_ticks

PLOTS = os.path.join(os.getcwd(), "plots")
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

ZOOM_COLOR = "tab:purple"  # tab:orange is taken by the 3D curve
ZOOM_LW = 1.0

_spec = importlib.util.spec_from_file_location(
    "stab3d2d", os.path.join(PLOTS, "3d2dStabilityCurveRe1000.py")
)
base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(base)


def plot(show_zoom, save_path, zoom_color="0.4", zoom_lw=0.5):
    frame = FigureFrame(
        frame_w=5.5, aspect_ratio=0.75, pad_l=4.0, pad_b=0.8, pad_t=0.11
    )
    fig, ax = frame.fig, frame.ax

    INS_W = 2.6
    INS_H = INS_W * frame.frame_h / frame.frame_w
    INS_LEFT = 0.55
    INS_BOTTOM = 0.8
    axins = frame.add_axes_cm(INS_LEFT, INS_BOTTOM, INS_W, INS_H)

    hopf2ds = average_curves(hopf2d_re_stable[1000], hopf_re_unstable[1000], npts=500)
    inst3ds = average_curves(
        instability3d_stable[1000], instability3d_unstable[1000], npts=10000
    )

    ax.set_xticks([0, 20, 40, 60, 80, 100, 120])
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_xlim([0, 130])
    ax.set_ylim([0, 4])
    axins.set_xlim([55, 85])
    axins.set_ylim([1.3, 1.75])

    zoom_kw = dict(fc="none", ec=zoom_color, lw=zoom_lw)
    zoom = TransformedBbox(axins.viewLim, ax.transData)
    artists = [ax.add_patch(BboxPatch(zoom, **zoom_kw))]
    for loc_ins, loc_zoom in ((1, 1), (4, 4)):
        connector = BboxConnector(
            axins.bbox, zoom, loc1=loc_ins, loc2=loc_zoom, zorder=5, **zoom_kw
        )
        connector.set_clip_on(False)
        artists.append(ax.add_patch(connector))
    ax.set_xlabel(r"$w$")
    ax.set_ylabel(r"$d$", rotation=0, labelpad=10)
    ax.grid(False)
    axins.grid(False)
    frame.jfm_ticks()
    jfm_ticks(axins)
    frame.axis_on_top()
    axis_on_top(axins)

    label_bg = dict(facecolor="white", edgecolor="none", pad=1.3)
    ax.yaxis.label.set_bbox(label_bg)
    for label in ax.get_yticklabels():
        label.set_bbox(label_bg)

    for a in [ax, axins]:
        base.makePatterns(a, inst3ds, hopf2ds)
        base.plotCurves(a, inst3ds, hopf2ds, exp_bypass_curve)
    base.addText(ax, axins)
    base.addLegend(ax)

    if not show_zoom:
        # hidden rather than removed, so the figure keeps exactly its layout
        axins.set_visible(False)
        for art in artists:
            art.set_visible(False)

    fig.savefig(save_path, format="pdf")
    print("Saved", os.path.normpath(save_path))


if __name__ == "__main__":
    plot(False, os.path.join(OUT_DIR, "stabCurveStep1.pdf"))
    plot(True, os.path.join(OUT_DIR, "stabCurveStep2.pdf"), ZOOM_COLOR, ZOOM_LW)
