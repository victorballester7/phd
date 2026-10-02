import os

import matplotlib
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from pp.colors import colors, mix_with_black
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import pandas as pd
from pp.figureFrame import FigureFrame

plt.style.use("plots/style/jfm.mplstyle")

MARKERS = ["o", "s", "^", "P", "*", "X", "h", "8"]
# 4 is the default marker size
MARKERSIZES = [4, 4, 4, 5, 6]  # tuned per glyph

runs = ["inc2dRe800", "inc2dRe1000", "inc2dRe1400", "inc2dRe2000", "inc2dRe3000"]
# runs = ["inc2dRe800"]

CMAP = matplotlib.colormaps.get_cmap("viridis")  # swap to "viridis" if preferred


def getRe(run: str):
    """
    extract the reynolds number from the runs
    """
    # Extract the part after "Re"
    re_str = run.split("Re")[-1]
    # Convert to integer
    re = int(re_str)
    return re


def plot_freq(df: pd.DataFrame, save_path: str):
    frame = FigureFrame(frame_w=5.6, aspect_ratio=0.85, pad_l=3.85, pad_b=0.8, pad_t=0.15)
    fig, ax = frame.fig, frame.ax

    ax.set_xlim([15, 95])
    ax.set_ylim([0.07, 0.28])
    x_tmp = np.linspace(10, 110, 100)
    ax.fill_between(
        x_tmp,
        0,
        0.127,
        color="maroon",
        alpha=0.2,
    )
    ax.plot(
        x_tmp,
        0.127 * np.ones_like(x_tmp),
        color="maroon",
        linestyle="--",
        linewidth=1.0,
    )
    ax.text(
        0.48,
        0.04,
        "Convectively unstable\nregion at " + r"${\mbox{\textit{Re}}}=1000$",
        transform=ax.transAxes,
        color="maroon",
        ha="left",
        va="bottom",
    )

    # Global depth range for a consistent colormap across all runs
    all_depths = np.unique(df.loc[df["run"].isin(runs), "d"])
    norm = mcolors.Normalize(vmin=1, vmax=4)

    for r, run in enumerate(runs):
        if run not in df["run"].values:
            continue

        re = getRe(run)
        df_run = df.loc[df["run"] == run]
        depths = np.unique(df_run["d"])
        m = MARKERS[r % len(MARKERS)]
        ms = MARKERSIZES[r % len(MARKERSIZES)]
        label_done = False

        for depth in depths:
            df_depth = df_run.loc[
                (df_run["d"] == depth) & (df_run["type"].isin(["limitcycle"]))
                # & (df_run["type"].isin(["equilibrium", "limitcycle"]))
            ]
            color = CMAP(norm(depth))

            # a visible edge: without it the dense cluster at small w merges
            color_edge = mix_with_black(color, alpha=0.45)

            ax.plot(
                df_depth["w"],
                df_depth["omega"],
                # df_depth["omega"],
                marker=m,
                markersize=ms,
                markerfacecolor=color,
                markeredgecolor=color_edge,
                markeredgewidth=0.5,
                color=color,
                linestyle="None",
                label=rf"$Re={re}$" if not label_done else "_nolegend_",
            )
            label_done = True

    # Colorbar for depth
    sm = cm.ScalarMappable(cmap=CMAP, norm=norm)
    sm.set_array([])

    cb = frame.colorbar(sm)
    cb.set_label(r"$d$", rotation=0, labelpad=5)

    # ax.set_xlim([10,110])
    ax.set_xlabel(r"$w$")
    ax.set_ylabel(r"$\omega$", rotation=0, labelpad=10)
    ax.set_xticks([20, 40, 60, 80])

    legend_handles = []

    for r, run in enumerate(runs):
        re = getRe(run)
        m = MARKERS[r % len(MARKERS)]
        ms = MARKERSIZES[r % len(MARKERSIZES)]

        legend_handles.append(
            Line2D(
                [0],
                [0],
                marker=m,
                markersize=ms,
                linestyle="None",
                # neutral grey: colour encodes the gap depth, so a coloured
                # legend marker would read as a particular depth
                markerfacecolor="0.65",
                markeredgecolor="0.25",
                markeredgewidth=0.5,
                label=rf"$Re={re}$",
            )
        )

    ax.grid(False)
    frame.axis_on_top()
    frame.jfm_ticks()

    # right of the plot like the other figures, but past the colorbar, its
    # tick labels and its label (CBAR_SPACE cm from the frame)
    CBAR_SPACE = 1.3
    ax.legend(
        handles=legend_handles,
        loc="center left",
        bbox_to_anchor=(1.05 + CBAR_SPACE / frame.frame_w, 0.5),
        handlelength=1.0,
        frameon=False,
    )

    fig.savefig(save_path, format="pdf")
    print(colors.OKGREEN + f"✓ Plot saved to: {save_path}" + colors.ENDC)


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_path, "../../../data/stability/stability.dat")
    save_path = os.path.join(script_path, "../../../images/freqUnstableMode.pdf")
    df = pd.read_csv(
        data_dir,
        sep=r"\s+",
        comment="#",
        names=["w", "d", "sigma", "omega", "type", "run"],
    )
    plot_freq(df, save_path)


if __name__ == "__main__":
    main()
