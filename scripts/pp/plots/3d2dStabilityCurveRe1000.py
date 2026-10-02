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
from mpl_toolkits.axes_grid1.inset_locator import BboxConnector, BboxPatch
from matplotlib.transforms import TransformedBbox
from matplotlib.backends import backend_pdf
from pp.figureFrame import FigureFrame, axis_on_top, jfm_ticks


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
# The five regions are drawn as vector hatched polygons bounded by the mean
# stability curves.  Regions II/III are split at the crossing of the two mean
# curves, and III/IV at X_SEP_III_IV.
#
# The PDF backend writes hatch patterns in opaque RGB (alpha is dropped), so
# the hatch colour is pre-blended with the white background instead.

# Single-character patterns, i.e. the sparse form of each: doubling a hatch
# character doubles its density, and at this figure size the doubled set read
# as one uniform texture instead of five distinguishable regions.
HATCH_PATTERNS = ["/", "o", "*", "O", "x"]
HATCH_COLOR = "tab:blue"
ALPHA_PATTERN = 0.48  # opacity of the hatch, emulated by blending with white
X_SEP_III_IV = 70.0  # w/d* at which region III gives way to region IV
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


def makePatterns(ax, inst3ds, hopf2ds):
    mpl.rcParams["hatch.linewidth"] = 0.5

    inst3d = inst3ds[2]
    hopf2d = hopf2ds[2]
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()

    # crossing of the two mean curves -> separates regions II and III
    xs = np.linspace(0.0, np.max(inst3d[:, 0]), 2000)
    gap = _curve_y(xs, inst3d) - _curve_y(xs, hopf2d)
    x_cross = xs[np.argmin(np.abs(gap))]

    def bounds(xa, xb):
        x = np.linspace(xa, xb, N_REGION_PTS)
        y2d = _curve_y(x, hopf2d)
        y3d = _curve_y(x, inst3d)
        return x, np.minimum(y2d, y3d), np.maximum(y2d, y3d)

    x, lo, hi = bounds(x0, x1)
    regions = [
        (x, np.full_like(x, y0), lo),  # I   below both
        bounds(x0, x_cross),  # II  between, left of the crossing
        bounds(x_cross, X_SEP_III_IV),  # III
        bounds(X_SEP_III_IV, x1),  # IV
        (x, hi, np.full_like(x, y1)),  # V   above both
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
        color="tab:brown",
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


def addText(ax, axins):
    # black, not the hatch colour: a blue numeral sitting on blue hatching
    # reads as part of the texture
    textcolor = "tab:blue"
    textcolor = mix_with_black(textcolor, 0.35)
    textweight = "bold"
    textfontsize = [14, 11]

    labels = ["I", "II", "III", "IV", "V"]

    positionsAx = {
        labels[0]: (0.25, 0.16),
        labels[1]: (0.14, 0.6),
        labels[2]: (0.48, 0.34),
        labels[3]: (0.8, 0.38),
        labels[4]: (0.535, 0.755),
    }
    positionsAxins = {
        labels[0]: (0.17, 0.13),
        labels[1]: (np.nan, np.nan),  # not visible in the inset
        labels[2]: (0.35, 0.4),
        labels[3]: (0.8, 0.45),
        labels[4]: (0.485, 0.82),
    }

    for a, p, t in zip([ax, axins], [positionsAx, positionsAxins], textfontsize):
        for label, (x, y) in p.items():
            a.text(
                x,
                y,
                label,
                transform=a.transAxes,
                fontsize=t,
                color=textcolor,
                fontweight=textweight,
                ha="center",
            )

    
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
            Line2D([0], [0], color="tab:brown", linestyle="--"),
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
        frame_w=5.5, aspect_ratio=0.75, pad_l=4.0, pad_b=0.8, pad_t=0.11
    )
    fig, ax = frame.fig, frame.ax

    # zoom panel in the left padding, vertically centred on the frame
    INS_W = 2.6
    INS_H = INS_W * frame.frame_h / frame.frame_w
    INS_LEFT = 0.55  # cm from the left figure edge (room for its tick labels)
    INS_BOTTOM = 0.8
    axins = frame.add_axes_cm(INS_LEFT, INS_BOTTOM, INS_W, INS_H)

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
    axins.set_xlim([55, 85])
    axins.set_ylim([1.3, 1.75])

    # zoomed region on the main axes, joined to the panel on its left:
    # mark_inset can only join matching corners, so the connectors are built
    # by hand (panel upper/lower right -> region upper/lower left)
    zoom_kw = dict(fc="none", ec="0.4", lw=0.5)
    zoom = TransformedBbox(axins.viewLim, ax.transData)
    ax.add_patch(BboxPatch(zoom, **zoom_kw))
    # the connectors belong to ax, below its axis artists (zorder 10 from
    # axis_on_top), so they pass under the y tick labels and the y label
    for loc_ins, loc_zoom in ((1, 1), (4, 4)):
        connector = BboxConnector(
            axins.bbox, zoom, loc1=loc_ins, loc2=loc_zoom, zorder=5, **zoom_kw
        )
        connector.set_clip_on(False)
        ax.add_patch(connector)
    ax.set_xlabel(r"$w$")
    ax.set_ylabel(r"$d$", rotation=0, labelpad=10)
    # remove grid
    ax.grid(False)
    axins.grid(False)
    frame.jfm_ticks()
    jfm_ticks(axins)
    # draw ticks and spines above the patterns, bands and curves
    frame.axis_on_top()
    axis_on_top(axins)

    # white behind the y labels, so the connectors read as passing beneath
    # them rather than through them
    label_bg = dict(facecolor="white", edgecolor="none", pad=1.3)
    ax.yaxis.label.set_bbox(label_bg)
    for label in ax.get_yticklabels():
        label.set_bbox(label_bg)

    for a in [ax, axins]:
        makePatterns(a, inst3ds, hopf2ds)
        plotCurves(a, inst3ds, hopf2ds, exp_bypass_curve)
    addText(ax, axins)
    addLegend(ax)

    # Save with fixed figure size (no bbox_inches="tight" to preserve exact dimensions)
    fig.savefig(save_path, format="pdf")
    print(colors.OKGREEN + f"✓ Plot saved to: {save_path}" + colors.ENDC)


if __name__ == "__main__":
    main()
