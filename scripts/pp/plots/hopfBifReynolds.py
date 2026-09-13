import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerBase
from matplotlib.patches import Rectangle

from pp.colors import colors
from pp.hopfDataPoints import hopf2d_re_stable, hopf_re_unstable
from pp.filterData import average_curves

plt.style.use("plots/style/jfm.mplstyle")

deltaStarsStartDomain = {
    800: 0.885946693021087,
    1000: 0.909901771786384,
    1400: 0.936529014995569,
    2000: 0.956012874991754,
    3000: 0.970896704821545,
}


# =============================================================================
# Geometry utilities
# =============================================================================
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


def plot_hopf_re(save_path):

    latex_width_cm = 8.0
    fig_width_in = latex_width_cm / 2.54

    _, ax = plt.subplots(
        figsize=(fig_width_in, fig_width_in * 0.75),
        constrained_layout=False,
    )

    # ax.set_position([
    #     0.17,  # left
    #     0.16,  # bottom
    #     0.66,  # width
    #     0.76,  # height
    # ])
    # Fix axis limits early so fill_between has correct bounds
    ax.set_xlim([0, 134])
    ax.set_ylim([0, 4.1])

    linestyles = [
        "-",  # solid
        (0, (1, 1)),  # very dense dots
        (0, (5, 5)),  # medium dashed
        (0, (10, 3)),  # long dash
        (0, (3, 1, 1, 1)),  # dash-dot-dot pattern
    ]

    legend_elements = []
    alpha = 0.3
    
    n = 5
    cmap_tmp = plt.get_cmap("brg_r")
    cmap = cmap_tmp(np.linspace(0.1, 0.9, n))

    if len(hopf2d_re_stable.keys()) != len(hopf_re_unstable.keys()):
        raise ValueError(
            "Stable and unstable branches must have the same number of Reynolds numbers."
        )

    for i, re in enumerate(hopf2d_re_stable.keys()):
        if re not in hopf_re_unstable:
            raise ValueError(f"Reynolds number {re} is missing in the unstable branch.")

        stable = hopf2d_re_stable[re]
        unstable = hopf_re_unstable[re]

        stable_smooth, unstable_smooth, middle_curve = average_curves(
            stable,
            unstable,
            npts=500,
        )

        # print(f"Re = {re}")
        # print(stable_smooth[:10], unstable_smooth[:10], middle_curve[:10])

        color = cmap[i]
        # color = plt.get_cmap("tab10")(i)
        ls = linestyles[i % len(linestyles)]

        ax.fill(
            np.r_[
                stable_smooth[:, 0],
                unstable_smooth[::-1, 0],
            ],
            np.r_[
                stable_smooth[:, 1],
                unstable_smooth[::-1, 1],
            ],
            color=color,
            alpha=alpha,
            linewidth=0,
            zorder=1,
        )

        # scaling_factor = deltaStarsStartDomain[re] # just to try if curves collapse each other
        scaling_factor = 1

        ax.plot(
            middle_curve[:, 0] * scaling_factor,
            middle_curve[:, 1] * scaling_factor,
            color=color,
            linestyle=ls,
            zorder=3,
        )
        line = Line2D([0], [0], color=color, linestyle=ls)

        patch = Rectangle(
            (0, 0),
            1,
            1,
            facecolor=color,
            alpha=alpha,
            edgecolor="none",
        )

        legend_elements.append((line, patch, rf"$\mathit{{Re}} = {re}$"))

    ax.set_xticks([0, 20, 40, 60, 80, 100, 120])

    ax.set_xlim(0, 130)
    ax.set_ylim(0, 4.1)

    ax.set_xlabel(r"$w/\delta^*$")
    ax.set_ylabel(r"$d/\delta^*$", rotation=0, labelpad=10)

    handles = [(line, patch) for line, patch, _ in legend_elements]
    labels = [label for _, _, label in legend_elements]

    ax.legend(
        handles,
        labels,
        # loc="center left",
        # bbox_to_anchor=(1.02, 0.5),
        handler_map={h: HandlerPatchLine() for h in handles},  # ← instance keys
        # borderaxespad=0,
        # frameon=True,
    )

    plt.savefig(
        save_path,
        format="pdf",
        bbox_inches="tight",
    )

    print(colors.OKGREEN + f"✓ Plot saved to: {save_path}" + colors.ENDC)


def main():

    script_path = os.path.dirname(os.path.abspath(__file__))

    output_path = os.path.join(
        script_path,
        "../../../images/hopfBifReynolds.pdf",
    )

    plot_hopf_re(output_path)


if __name__ == "__main__":
    main()
