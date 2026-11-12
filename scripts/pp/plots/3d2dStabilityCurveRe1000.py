import numpy as np
import matplotlib.pyplot as plt
import os
from pp.colors import colors, mix_with_black
from pp.hopfDataPoints import (
    hopf2d_re_stable,
    hopf_re_unstable,
    instability3d_stable,
    instability3d_unstable,
    exp_bypass_curve,
)
from pp.filterData import average_curves
from matplotlib.legend_handler import HandlerBase
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import matplotlib as mpl
from matplotlib.backends import backend_pdf
from pp.figureFrame import FigureFrame


plt.style.use("plots/style/jfm.mplstyle")


class HandlerPatchLine(HandlerBase):
    def create_artists(
        self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans
    ):
        line, patch = orig_handle

        p = Rectangle(
            (xdescent, ydescent),
            width,
            height,
            facecolor=patch.get_facecolor(),
            edgecolor="none",
            alpha=patch.get_alpha(),
            transform=trans,
        )

        l = Line2D(
            [xdescent, xdescent + width],
            [ydescent + height / 2, ydescent + height / 2],
            transform=trans,
        )
        l.update_from(line)
        l.set_data(
            [xdescent, xdescent + width],
            [ydescent + height / 2, ydescent + height / 2],
        )
        l.set_transform(trans)
        l.set_solid_capstyle("butt")  # ← prevents cap bleed beyond endpoints

        return [p, l]


# --- Region patterns -----------------------------------------------------
#
# The four regions are drawn as vector hatched polygons bounded by the mean
# stability curves.  Regions II/III are split at the crossing of the two mean
# curves.
#
# The PDF backend writes hatch patterns in opaque RGB (alpha is dropped), so
# the hatch colour is pre-blended with the white background instead.

# Single-character patterns, i.e. the sparse form of each: doubling a hatch
# character doubles its density, and at this figure size the doubled set read
# as one uniform texture instead of four distinguishable regions.
HATCH_PATTERNS = ["/", "o", "*", "x"]
HATCH_COLOR = "tab:blue"
ALPHA_PATTERN = 0.48  # opacity of the hatch, emulated by blending with white
N_REGION_PTS = 2000  # resolution of the region outlines


def _writeHatches(self):
    """PdfFile.writeHatches, but filling closed hatch shapes (e.g. the "*"
    stars) with the hatch colour: upstream leaves the fill colour unset when
    the patch has no facecolor, so they come out black."""
    hatchDict = dict()
    sidelen = 72.0
    Op = backend_pdf.Op
    for hatch_style, name in self._hatch_patterns.items():
        ob = self.reserveObject("hatch pattern")
        hatchDict[name] = ob
        res = {
            "Procsets": [
                backend_pdf.Name(x) for x in "PDF Text ImageB ImageC ImageI".split()
            ]
        }
        self.beginStream(
            ob.id,
            None,
            {
                "Type": backend_pdf.Name("Pattern"),
                "PatternType": 1,
                "PaintType": 1,
                "TilingType": 1,
                "BBox": [0, 0, sidelen, sidelen],
                "XStep": sidelen,
                "YStep": sidelen,
                "Resources": res,
                "Matrix": [1, 0, 0, 1, 0, self.height * 72],
            },
        )
        stroke_rgb, fill_rgb, hatch, lw = hatch_style
        self.output(*stroke_rgb[:3], Op.setrgb_stroke)
        if fill_rgb is not None:
            self.output(
                *fill_rgb[:3],
                Op.setrgb_nonstroke,
                0,
                0,
                sidelen,
                sidelen,
                Op.rectangle,
                Op.fill,
            )
        self.output(*stroke_rgb[:3], Op.setrgb_nonstroke)  # <- the fix
        self.output(lw, Op.setlinewidth)
        self.output(
            *self.pathOperations(
                backend_pdf.Path.hatch(hatch),
                backend_pdf.Affine2D().scale(sidelen),
                simplify=False,
            )
        )
        self.output(Op.fill_stroke)
        self.endStream()
    self.writeObject(self.hatchObject, hatchDict)


backend_pdf.PdfFile.writeHatches = _writeHatches


def _curve_y(x, curve):
    """y(x) of an (N, 2) curve, tolerant of the tiny non-monotonicity left by
    the Savitzky-Golay smoothing in average_curves."""
    order = np.argsort(curve[:, 0])
    return np.interp(x, curve[order, 0], curve[order, 1])


def _curve_crossing(inst3d, hopf2d):
    """w at which the two mean stability curves cross (regions II | III)."""
    xs = np.linspace(0.0, np.max(inst3d[:, 0]), 2000)
    gap = _curve_y(xs, inst3d) - _curve_y(xs, hopf2d)
    return xs[np.argmin(np.abs(gap))]


def makePatterns(ax, inst3ds, hopf2ds):
    mpl.rcParams["hatch.linewidth"] = 0.5

    inst3d = inst3ds[2]
    hopf2d = hopf2ds[2]
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()

    # crossing of the two mean curves -> separates regions II and III
    x_cross = _curve_crossing(inst3d, hopf2d)

    def bounds(xa, xb):
        x = np.linspace(xa, xb, N_REGION_PTS)
        y2d = _curve_y(x, hopf2d)
        y3d = _curve_y(x, inst3d)
        return x, np.minimum(y2d, y3d), np.maximum(y2d, y3d)

    x, lo, hi = bounds(x0, x1)
    regions = [
        (x, np.full_like(x, y0), lo),  # I   below both
        bounds(x0, x_cross),  # II  between, left of the crossing
        bounds(x_cross, x1),  # III between, right of the crossing
        (x, hi, np.full_like(x, y1)),  # IV  above both
    ]

    rgb = np.asarray(mpl.colors.to_rgb(HATCH_COLOR))
    edgecolor = ALPHA_PATTERN * rgb + (1.0 - ALPHA_PATTERN)  # over white
    for (x, ya, yb), hatch_pattern in zip(regions, HATCH_PATTERNS):
        ax.fill_between(
            x,
            ya,
            yb,
            facecolor="none",
            edgecolor=edgecolor,
            linewidth=0,
            hatch=hatch_pattern,
            zorder=1,
        )


def plotCurves(ax, inst3ds, hopf2ds, exp_bypass_curve):
    inst3d_st, inst3d_un, inst3d = inst3ds
    hopf2d_st, hopf2d_un, hopf2d = hopf2ds

    ax.plot(
        hopf2d[:, 0],
        hopf2d[:, 1],
        linestyle="-",
        color="tab:green",
        zorder=2,
    )

    ax.plot(
        inst3d[:, 0],
        inst3d[:, 1],
        linestyle="-.",
        color="tab:orange",
        zorder=3,
    )

    ax.plot(
        exp_bypass_curve[:, 0],
        exp_bypass_curve[:, 1],
        linestyle="--",
        color="k",
        zorder=3,
    )

    ax.fill(
        np.r_[hopf2d_st[:, 0], hopf2d_un[::-1, 0]],
        np.r_[hopf2d_st[:, 1], hopf2d_un[::-1, 1]],
        color="tab:green",
        alpha=0.2,
        linewidth=0,
        zorder=1,
    )

    ax.fill(
        np.r_[inst3d_st[:, 0], inst3d_un[::-1, 0]],
        np.r_[inst3d_st[:, 1], inst3d_un[::-1, 1]],
        color="tab:orange",
        alpha=0.2,
        linewidth=0,
        zorder=1,
    )


def addText(ax, inst3d, hopf2d):
    # black, not the hatch colour: a blue numeral sitting on blue hatching
    # reads as part of the texture
    textcolor = "tab:blue"
    textcolor = mix_with_black(textcolor, 0.35)
    kw = dict(fontsize=14, color=textcolor, fontweight="bold", ha="center")

    for label, (x, y) in {
        "I": (0.25, 0.16),
        "II": (0.14, 0.6),
        "IV": (0.535, 0.755),
    }.items():
        ax.text(x, y, label, transform=ax.transAxes, **kw)

    # III: centre of the thin strip between the two curves, right of their
    # crossing, in data coordinates
    x_cross = _curve_crossing(inst3d, hopf2d)
    x_mid = 0.5 * (x_cross + ax.get_xlim()[1])
    y_mid = 0.5 * (_curve_y(x_mid, inst3d) + _curve_y(x_mid, hopf2d))
    ax.text(x_mid, y_mid, "III", va="center_baseline", **kw)

    
    ax.text(
        0.8,
        0.3,
        r"$d_{\mathrm{2D}}$",
        transform=ax.transAxes,
        color="tab:green",
        ha="center",
    )
    ax.text(
        0.24,
        0.36,
        r"$d_{\mathrm{3D}}$",
        transform=ax.transAxes,
        color="tab:orange",
        ha="center",
    )


def addLegend(ax):
    legend_elements = [
        (
            Line2D([0], [0], color="tab:green", linestyle="-"),
            Rectangle((0, 0), 1, 1, facecolor="tab:green", alpha=0.2),
            "2D linear stability\nboundary",
        ),
        (
            Line2D([0], [0], color="tab:orange", linestyle="-."),
            Rectangle((0, 0), 1, 1, facecolor="tab:orange", alpha=0.2),
            "3D linear stability\nboundary",
        ),
        (
            Line2D([0], [0], color="k", linestyle="--"),
            Rectangle((0, 0), 1, 1, facecolor="white", alpha=0.0),
            "Experimental bypass\ntransition boundary\n(Crouch $\\it{et\\ al.}$ 2022)",
        ),
    ]

    handles = [(line, patch) for line, patch, _ in legend_elements]
    labels = [label for _, _, label in legend_elements]

    ax.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        handler_map={h: HandlerPatchLine() for h in handles},  # ← instance keys
        frameon=False,
    )


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(
        script_path, "../../../images/3d2dStabilityCurveRe1000.pdf"
    )

    frame = FigureFrame(
        frame_w=5.5, aspect_ratio=0.75, pad_l=4.0, pad_b=0.8, pad_t=0.11  # pad_r = pad_l: frame centred, legend in pad_r
    )
    fig, ax = frame.fig, frame.ax

    hopf2ds = average_curves(hopf2d_re_stable[1000], hopf_re_unstable[1000], npts=500)

    inst3ds = average_curves(
        instability3d_stable[1000], instability3d_unstable[1000], npts=10000
    )

    # NOTE: set_ticks widens the view interval to span the ticks, so the
    # limits must be set *after* the ticks or they get silently overridden.
    ax.set_xticks([0, 20, 40, 60, 80, 100, 120])
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_xlim([0, 130])
    ax.set_ylim([0, 4])

    ax.set_xlabel(r"$w$")
    ax.set_ylabel(r"$d$", rotation=0, labelpad=10)
    ax.grid(False)
    frame.jfm_ticks()
    # draw ticks and spines above the patterns, bands and curves
    frame.axis_on_top()

    makePatterns(ax, inst3ds, hopf2ds)
    plotCurves(ax, inst3ds, hopf2ds, exp_bypass_curve)
    addText(ax, inst3ds[2], hopf2ds[2])
    addLegend(ax)

    # Save with fixed figure size (no bbox_inches="tight" to preserve exact dimensions)
    fig.savefig(save_path, format="pdf")
    print(colors.OKGREEN + f"✓ Plot saved to: {save_path}" + colors.ENDC)


if __name__ == "__main__":
    main()
