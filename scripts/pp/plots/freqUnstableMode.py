import os

import matplotlib
from matplotlib.lines import Line2D
from pp.colors import colors, mix_with_black
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import pandas as pd

plt.style.use("plots/style/jfm.mplstyle")

MARKERS = ["o", "s", "^", "P", "*", "X", "h", "8"]
# 4 is the default marker size
MARKERSIZES = [4, 4, 4, 5, 6]  # tuned per glyph

runs = ["inc2dRe800", "inc2dRe1000", "inc2dRe1400", "inc2dRe2000", "inc2dRe3000"]
# runs = ["inc2dRe800"]

CMAP = matplotlib.colormaps.get_cmap("viridis")  # swap to "viridis" if preferred
# CMAP = matplotlib.colormaps.get_cmap("Blues")  # swap to "viridis" if preferred


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
    latex_width_cm = 9.0
    fig_width_in = latex_width_cm / 2.54

    fig, ax = plt.subplots(
        figsize=(fig_width_in, fig_width_in * 0.75),
        constrained_layout=False,
    )

    ax.set_xlim([15, 95])
    ax.set_ylim([0.07, 0.28])
    x_tmp = np.linspace(10, 110, 100)
    ax.fill_between(
        x_tmp,
        0,
        0.127,
        color="tab:red",
        alpha=0.2,
    )
    ax.plot(
        x_tmp,
        0.127 * np.ones_like(x_tmp),
        color="darkred",
        linestyle="--",
        linewidth=1.0,
    )
    ax.text(
        0.48,
        0.015,
        "Convectively unstable\nregion at " + r"${\mbox{\textit{Re}}}=1000$",
        transform=ax.transAxes,
        color="darkred",
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

            color_edge = mix_with_black(color, alpha=0.1)

            ax.plot(
                df_depth["w"],
                df_depth["omega"],
                # df_depth["omega"],
                marker=m,
                markersize=ms,
                markerfacecolor=color,
                markeredgecolor=color_edge,
                color=color,
                linestyle="None",
                label=rf"$Re={re}$" if not label_done else "_nolegend_",
            )
            label_done = True

    # Colorbar for depth
    sm = cm.ScalarMappable(cmap=CMAP, norm=norm)
    sm.set_array([])

    cb = fig.colorbar(sm, ax=ax)
    cb.set_label(r"$d/\delta^*$", rotation=0, labelpad=15)

    # ax.set_xlim([10,110])
    ax.set_xlabel(r"$w/\delta^*$")
    ax.set_ylabel(r"$\omega/(\frac{\delta^*}{u_\infty})$", rotation=0, labelpad=18)
    ax.set_xticks([20, 40, 60, 80])
    ax.legend()

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
                markerfacecolor=CMAP(norm(3.0)),  # dark blue
                markeredgecolor=mix_with_black(CMAP(norm(3.0)), alpha=0.1),
                label=rf"$Re={re}$",
            )
        )

    ax.legend(handles=legend_handles)

    plt.savefig(save_path, format="pdf", bbox_inches="tight")
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
