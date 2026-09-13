import numpy as np
import matplotlib.pyplot as plt
import os
from pp.colors import colors
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


def create_wiggle_curve(x0, y0, x_end, slope, amplitude, frequency, num_points=500):
    x_wiggle = np.linspace(x0, x_end, num_points)
    s_black = x_wiggle - x0
    y_wiggle = (
        y0
        + slope * s_black
        + amplitude * np.sin(frequency * np.pi * s_black / max(s_black[-1], 1.0))
    )
    return x_wiggle, y_wiggle


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(
        script_path, "../../../images/3d2dStabilityCurveRe1000.pdf"
    )

    latex_width_cm = 8.0
    fig_width_in = latex_width_cm / 2.54

    _, ax = plt.subplots(
        figsize=(fig_width_in, fig_width_in * 0.75),
        constrained_layout=False,
    )

    hopf2d_st, hopf2d_un, hopf2d = average_curves(
        hopf2d_re_stable[1000], hopf_re_unstable[1000], npts=500
    )

    inst3d_st, inst3d_un, inst3d = average_curves(
        instability3d_stable[1000], instability3d_unstable[1000], npts=10000
    )

    # Fix axis limits early so fill_between has correct bounds
    ax.set_xlim([0, 130])
    ax.set_ylim([0, 4])

    # filter of hopf2d, hopf2d_st and hopf2d_un to be less than 4 for the depth
    hopf2d_st = hopf2d_st[hopf2d_st[:, 1] < 4]
    hopf2d_un = hopf2d_un[hopf2d_un[:, 1] < 4]
    hopf2d = hopf2d[hopf2d[:, 1] < 4]

    #### curves with filling

    # wiggle
    x01 = 48
    x_wiggle1, y_wiggle1 = create_wiggle_curve(
        x0=x01,
        y0=1.7,
        x_end=130,
        slope=0.005,
        amplitude=0.05,
        frequency=5,
        num_points=500,
    )

    x02 = 72
    x_end2 = 69.8
    x_wiggle2, y_wiggle2 = create_wiggle_curve(
        x0=x02,
        y0=1.5,
        x_end=x_end2,
        slope=-0.1,
        amplitude=0.05,
        frequency=0.6,
        num_points=100,
    )

    # Plot the black uncertainty boundary.
    for x, y in zip([x_wiggle1, x_wiggle2], [y_wiggle1, y_wiggle2]):
        ax.plot(
            x,
            y,
            color="black",
            linewidth=0,
            zorder=6,
        )

    # filter values of inst3d_st and inst3d_un to be less than width = 50
    inst3d_st_filt = inst3d_st[inst3d_st[:, 0] < 50]
    inst3d_un_filt = inst3d_un[inst3d_un[:, 0] < 50]

    w_regionI = np.concatenate([[0], inst3d_st[:, 0]])
    d_2_regionI = np.concatenate([[4], inst3d_st[:, 1]])
    d_1_regionII = np.interp(inst3d_un_filt[:, 0], hopf2d_st[:, 0], hopf2d_st[:, 1])
    d_2_regionII = np.interp(inst3d_un_filt[:, 0], inst3d_un[:, 0], inst3d_un[:, 1])

    x_regionIV = np.where(x_wiggle1 >= x_end2)[0]
    y_1_regionIV = y_wiggle1[x_regionIV]
    x_regionIV = x_wiggle1[x_regionIV]
    y_2_regionIV_1 = np.interp(
        x_regionIV[np.where(x_regionIV <= x02)[0]], x_wiggle2[::-1], y_wiggle2[::-1]
    )
    y_2_regionIV_2 = np.interp(
        x_regionIV[np.where(x_regionIV > x02)[0]], hopf2d_un[:, 0], hopf2d_un[:, 1]
    )
    y_2_regionIV = np.concatenate([y_2_regionIV_1, y_2_regionIV_2])

    x_regionV_1 = np.where(x_wiggle1 <= x_end2)[0]
    y_regionV_1 = y_wiggle1[x_regionV_1]
    x_regionV_1 = x_wiggle1[x_regionV_1]
    x_regionV_2 = x_wiggle2[::-1]
    y_regionV_2 = y_wiggle2[::-1]

    x_regionV = np.concatenate([x_regionV_1, x_regionV_2])
    y_1_regionV = np.concatenate([y_regionV_1, y_regionV_2])
    y_2_regionV = np.interp(x_regionV, hopf2d_un[:, 0], hopf2d_un[:, 1])
    hopf2d_un_filt = hopf2d_un[hopf2d_un[:, 0] < x01]
    d_2_regionIII = np.concatenate([hopf2d_un_filt[:, 1], y_wiggle1])
    w_region_top = np.concatenate([hopf2d_un_filt[:, 0], x_wiggle1])
    d_2_regionIII = np.interp(hopf2d_un[:, 0], w_region_top, d_2_regionIII)
    #### filling of different regions
    hatch_patterns = ["//", "O.", "xx", "oo", "**"]
    alpha_fb = 0.3

    mpl.rcParams["hatch.linewidth"] = 0.7

    ax.fill_between(
        w_regionI,
        0,
        d_2_regionI,
        facecolor="white",
        edgecolor="tab:blue",
        linewidth=0,
        alpha=alpha_fb,
        hatch=hatch_patterns[0],
        zorder=1,
    )

    ax.fill_between(
        inst3d_un_filt[:, 0],
        d_1_regionII,
        d_2_regionII,
        facecolor="white",
        edgecolor="tab:blue",
        linewidth=0,
        alpha=alpha_fb,
        hatch=hatch_patterns[1],
        zorder=1,
    )

    ax.fill_between(
        hopf2d_un[:, 0],
        4,
        d_2_regionIII,
        facecolor="white",
        edgecolor="tab:blue",
        linewidth=0.5,
        alpha=alpha_fb,
        hatch=hatch_patterns[2],
        zorder=1,
    )

    ax.fill_between(
        x_regionIV,
        y_1_regionIV,
        y_2_regionIV,
        facecolor="white",
        edgecolor="tab:blue",
        linewidth=0,
        alpha=alpha_fb,
        hatch=hatch_patterns[3],
        zorder=1,
    )
    # print(x_regionIV, y_1_regionIV, y_2_regionIV)

    ax.fill_between(
        x_regionV,
        y_1_regionV,
        y_2_regionV,
        facecolor="white",
        edgecolor="tab:blue",
        linewidth=0,
        alpha=alpha_fb,
        hatch=hatch_patterns[4],
        zorder=1,
    )
    ax.plot(
        x_wiggle2,
        y_wiggle2,
        color="tab:blue",
        linewidth=0.5,
        alpha=alpha_fb,
        zorder=1,
    )

    ax.plot(
        hopf2d[:, 0],
        hopf2d[:, 1],
        linestyle="-",
        label="2D linear stability\nboundary",
        color="tab:green",
        zorder=2,
    )


    ax.plot(
        inst3d[:, 0],
        inst3d[:, 1],
        linestyle="-.",
        label="3D linear stability\nboundary",
        color="tab:orange",
        zorder=3,
    )

    ax.plot(
        exp_bypass_curve[:, 0],
        exp_bypass_curve[:, 1],
        linestyle="--",
        color="tab:brown",
        label="Experimental bypass\ntransition boundary\n(Crouch $\\it{et\\ al.}$ 2022)",
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
        np.r_[inst3d_st_filt[:, 0], inst3d_un_filt[::-1, 0]],
        np.r_[inst3d_st_filt[:, 1], inst3d_un_filt[::-1, 1]],
        color="tab:orange",
        alpha=0.2,
        linewidth=0,
        zorder=1,
    )

    ax.set_xticks([0, 20, 40, 60, 80, 100, 120])
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_xlabel(r"$w/\delta^*$")
    ax.set_ylabel(r"$d/\delta^*$", rotation=0, labelpad=10)

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

    textcolor = "tab:blue"
    textfontsize = 12

    ax.text(
        0.25,
        0.15,
        "I",
        transform=ax.transAxes,
        fontsize=textfontsize,
        color=textcolor,
        ha="center",
    )
    ax.text(
        0.15,
        0.6,
        "II",
        transform=ax.transAxes,
        fontsize=textfontsize,
        color=textcolor,
        ha="center",
    )
    ax.text(
        0.23,
        0.85,
        "V",
        transform=ax.transAxes,
        fontsize=textfontsize,
        color=textcolor,
        ha="center",
    )
    ax.text(
        0.7,
        0.4,
        "IV",
        transform=ax.transAxes,
        fontsize=textfontsize,
        color=textcolor,
        ha="center",
    )
    ax.text(
        0.51,
        0.36,
        "III",
        transform=ax.transAxes,
        fontsize=textfontsize,
        color=textcolor,
        ha="center",
    )

    # ax.plot([14, 24], [4, 2], color="tab:purple", marker="*", linestyle="none", markersize=8)
    # ax.plot([19, 40], [4, 2], color="tab:cyan", marker="*", linestyle="none", markersize=8)
    # ax.plot([80], [1.5], color="tab:red", marker="*", linestyle="none", markersize=8)

    handles = [(line, patch) for line, patch, _ in legend_elements]
    labels = [label for _, _, label in legend_elements]

    ax.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        handler_map={h: HandlerPatchLine() for h in handles},  # ← instance keys
    )

    # Save with fixed figure size (no bbox_inches="tight" to preserve exact dimensions)
    plt.savefig(save_path, format="pdf")
    print(colors.OKGREEN + f"✓ Plot saved to: {save_path}" + colors.ENDC)


if __name__ == "__main__":
    main()
