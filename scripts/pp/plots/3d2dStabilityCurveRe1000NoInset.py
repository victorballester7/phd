"""3d2dStabilityCurveRe1000 without the zoom inset, so the map fills the figure
(used e.g. as a thumbnail in slide outlines)."""

import os
import importlib.util
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

# the module name starts with a digit, so it cannot be imported normally
_spec = importlib.util.spec_from_file_location(
    "stab3d2d",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "3d2dStabilityCurveRe1000.py"),
)
base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(base)


def addText(ax):
    textcolor = mix_with_black("tab:blue", 0.35)
    positions = {
        "I": (0.25, 0.16),
        "II": (0.14, 0.6),
        "III": (0.48, 0.34),
        "IV": (0.8, 0.38),
        "V": (0.535, 0.755),
    }
    for label, (x, y) in positions.items():
        ax.text(
            x,
            y,
            label,
            transform=ax.transAxes,
            fontsize=14,
            color=textcolor,
            fontweight="bold",
            ha="center",
        )
    ax.text(0.8, 0.3, r"$d_{\mathrm{2D}}$", transform=ax.transAxes, color="tab:green", ha="center")
    ax.text(0.24, 0.36, r"$d_{\mathrm{3D}}$", transform=ax.transAxes, color="tab:orange", ha="center")


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(
        script_path, "../../../images/3d2dStabilityCurveRe1000NoInset.pdf"
    )

    frame = FigureFrame(
        frame_w=5.5, aspect_ratio=0.75, pad_l=0.9, pad_r=4.0, pad_b=0.8, pad_t=0.11
    )
    fig, ax = frame.fig, frame.ax

    hopf2ds = average_curves(hopf2d_re_stable[1000], hopf_re_unstable[1000], npts=500)
    inst3ds = average_curves(
        instability3d_stable[1000], instability3d_unstable[1000], npts=10000
    )

    # set_ticks widens the view interval, so set the limits afterwards
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
    base.plotCurves(ax, inst3ds, hopf2ds, exp_bypass_curve)
    addText(ax)
    base.addLegend(ax)

    fig.savefig(save_path, format="pdf")
    print(colors.OKGREEN + f"✓ Plot saved to: {save_path}" + colors.ENDC)


if __name__ == "__main__":
    main()
