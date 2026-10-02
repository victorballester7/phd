import os
import numpy as np
import matplotlib.pyplot as plt
from pp.colors import mix_with_black, colors
from pp.figureFrame import FigureFrame
import pandas as pd
from pp.filterData import average_curves
from pp.hopfDataPoints import hopf2d_re_stable, hopf_re_unstable, exp_bypass_curve


# plt.style.use('plots/style/tsfp.mplstyle')
plt.style.use("plots/style/jfm.mplstyle")


MARKERS = ["o", "^", "*"]
MARKERSIZE = [
    plt.rcParams["lines.markersize"],
    plt.rcParams["lines.markersize"],
    plt.rcParams["lines.markersize"] * 1.2,
]
MARKER_COLORS = ["tab:blue", "gold", "firebrick"]

# labels as dictionary
LABELS = {
    "equilibrium": "Equilibrium",
    "limitcycle": "Limit cycle",
    "chaoticattractor": "Chaotic state",
}

def plot_stability_diagram(df, save_path):
    frame = FigureFrame(frame_w=5.8, aspect_ratio=0.75, pad_l=3.85, pad_b=0.8,pad_t=0.05)
    fig, ax = frame.fig, frame.ax

    # Fix axis limits early so fill_between has correct bounds
    ax.set_xlim([0, 134])
    ax.set_ylim([0, 4.1])

    inc2dData = df.loc[df["run"] == "inc2dRe1000"].to_numpy()

    for m, ms, c, lab in zip(MARKERS, MARKERSIZE, MARKER_COLORS, LABELS.keys()):
        subset = inc2dData[inc2dData[:, 4] == lab]
        edge_color = mix_with_black(c, alpha=0.3)
        ax.plot(
            subset[:, 0],
            subset[:, 1],
            marker=m,
            markersize=ms,
            markerfacecolor=c,
            markeredgecolor=edge_color,
            alpha=0.5,
            linestyle="none",
            label=LABELS[lab],
            zorder=3,
        )

    _, _, hopf2d = average_curves(
        hopf2d_re_stable[1000], hopf_re_unstable[1000], npts=500
    )

    # Sort by w for fill_between
    w_hopf = hopf2d[:, 0]
    d_hopf = hopf2d[:, 1]

    unstable_color = (0.85, 0.30, 0.05)  # orange-red; tweak as needed

    # Extend the Hopf curve to the left and right axes edges
    # by prepending/appending boundary points
    w_plot = np.concatenate([w_hopf, [133]])
    d_plot = np.concatenate([d_hopf, [d_hopf[-1]]])  # Extend d to the last value
    w_extended = np.concatenate([[0], w_hopf, [140]])
    d_extended = np.concatenate(
        [[5], d_hopf, [d_hopf[-1]]]
    )  


    alpha_fb = 0.2

    # Stable: single fill_between, no fill_betweenx needed
    ax.fill_between(
        w_extended,
        0,
        d_extended,
        color="tab:blue",
        alpha=alpha_fb,
        zorder=1,
        label="Stable region",
    )

    # Unstable: same idea
    ax.fill_between(
        w_extended,
        d_extended,
        4.1,
        color=unstable_color,
        alpha=alpha_fb,
        zorder=1,
        label="Unstable region",
    )

    ax.plot(
        w_plot,
        d_plot,
        linestyle="-",
        color="tab:green",
        label="Linear stability\nboundary",
        zorder=4,
    )

    ax.plot(
        exp_bypass_curve[:, 0],
        exp_bypass_curve[:, 1],
        linestyle="--",
        color="tab:brown",
        label="Experimental bypass\ntransition boundary\n(Crouch $\\it{et\\ al.}$ 2022)",
        zorder=5,
    )

    ax.set_xticks([0, 20, 40, 60, 80, 100, 120])
    ax.set_xlabel(r"$w$")
    ax.set_ylabel(r"$d$", rotation=0, labelpad=10)

    # add location to the right of the plot for the legend
    # ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    ax.grid(False)

    # draw ticks and spines above the markers, filled regions and curves
    frame.axis_on_top()
    frame.jfm_ticks()
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),  # just outside the right edge of the axes
        frameon=False,
    )

    fig.savefig(
        save_path,
        format="pdf",
    )  

    print(colors.OKGREEN + f"✓ Plot saved to: {save_path}" + colors.ENDC)


def main():
    """Main execution function."""
    # Set up paths
    script_path = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_path, "../../../data/stability/stability.dat")
    output_path = os.path.join(
        script_path, "../../../images/inc2dStabilityCurveRe1000.pdf"
    )

    # Load data
    data_frames = pd.read_csv(
        data_dir,
        sep=r"\s+",
        comment="#",
        names=["w", "d", "sigma", "omega", "type", "run"],
    )

    # Create plot
    plot_stability_diagram(data_frames, output_path)


if __name__ == "__main__":
    main()
