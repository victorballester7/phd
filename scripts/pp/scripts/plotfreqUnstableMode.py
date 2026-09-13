import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import pandas as pd

MARKERS = ["o", "s", "^", "P", "*", "X", "h", "8"]
MARKERSIZES = [6, 6, 7, 8, 10, 8, 7, 8]  # tuned per glyph

runs = ["inc2dRe800", "inc2dRe1000", "inc2dRe1400", "inc2dRe2000", "inc2dRe3000"]
# runs = ["inc2dRe1000", "inc2dRe1400", "inc2dRe2000", "inc2dRe3000"]

CMAP = cm.get_cmap("inferno")  # swap to "viridis" if preferred


def plot_freq(df: pd.DataFrame):
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    # Global depth range for a consistent colormap across all runs
    all_depths = np.unique(df.loc[df["run"].isin(runs), "d"])
    d_min, d_max = all_depths.min(), all_depths.max()
    norm = mcolors.Normalize(vmin=d_min, vmax=d_max)

    for r, run in enumerate(runs):
        if run not in df["run"].values:
            continue

        df_run = df.loc[df["run"] == run]
        depths = np.unique(df_run["d"])
        marker = MARKERS[r % len(MARKERS)]
        ms = MARKERSIZES[r % len(MARKERSIZES)]
        label_done = False

        for depth in depths:
            df_depth = df_run.loc[
                (df_run["d"] == depth)
                & (df_run["type"].isin(["limitcycle"]))
                # & (df_run["type"].isin(["equilibrium", "limitcycle"]))
            ]
            color = CMAP(norm(depth))

            label = run if not label_done else "_nolegend_"
            label_done = True

            ax.plot(
                df_depth["w"],
                df_depth["omega"],
                marker=marker,
                markersize=ms,
                color=color,
                linestyle="None",
                label=label,
            )

    # Colorbar for depth
    sm = cm.ScalarMappable(cmap=CMAP, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.02)
    cb.set_label(r"$d/\delta^*$")

    ax.set_xlabel(r"$w/\delta^*$")
    ax.set_ylabel(r"$\omega$")
    ax.set_title("Frequency vs Width")
    ax.legend(title="Run")
    ax.grid()
    plt.show()


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_path, "../../../data/stability/stability.dat")
    df = pd.read_csv(
        data_dir,
        sep=r"\s+",
        comment="#",
        names=["w", "d", "sigma", "omega", "type", "run"],
    )
    plot_freq(df)


if __name__ == "__main__":
    main()
