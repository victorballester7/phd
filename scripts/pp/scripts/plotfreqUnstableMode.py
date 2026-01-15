import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

MARKERS = ["o", "s", "^", "P", "*", "X", "h", "8"]

runs = ["inc2d", "Ma0.2_2d", "Ma0.4_2d", "Ma0.6_2d", "Ma0.8_2d"]

depth_same_color = True

def plot_freq(df: pd.DataFrame):
    """
    plot the frequency vs width with different markers for df_1 and df_2,
    but show only ONE legend entry per depth.
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    print("Unique runs in data:", np.unique(df["run"]))

    for r, run in enumerate(np.unique(df["run"])):
        if run not in runs:
            continue

        print(f"Processing run: {run}")

        df_run = df.loc[df["run"] == run]

        depths = np.unique(df_run["d"])[::-1]

        label_done = False
        for i, depth in enumerate(depths):
            df_depth = df_run.loc[df_run["d"] == depth]


            # fix color, ith in the mpl10 color cycle
            if depth_same_color:
                color = plt.rcParams["axes.prop_cycle"].by_key()["color"][r % 10]
            else:
                color = plt.rcParams["axes.prop_cycle"].by_key()["color"][i % 10]

            # fix marker style based on run
            marker = MARKERS[r % len(MARKERS)]

            label = ""
            if not depth_same_color:
                label = f"Depth = {depth} (run: {run})"
            elif not label_done:
                label = f"Run: {run}"
                label_done = True
            else:
                label = "_nolegend_"

            ax.plot(
                df_depth["w"],
                df_depth["omega"],
                marker=marker,
                markersize=6,
                color=color,
                linestyle="None",
                # linestyle=linestyles[r % len(linestyles)],
                # Only label ONCE per depth (for df1)
                label=label,
            )

            ax.plot(
                df_depth["w"],
                df_depth["omega"],
                color=color,
                linewidth=1.0,
                marker=None,
                linestyle="-",
                label="_nolegend_",
            )

    ax.set_xlabel("w")
    ax.set_ylabel("ω")
    ax.set_title("Frequency vs Width")
    ax.legend()
    ax.grid()
    plt.show()


def main():
    """Main execution function."""
    # Set up paths
    script_path = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_path, "../../../data/stability/stability.dat")

    # Load data
    df = pd.read_csv(
        data_dir,
        sep=" ",
        comment="#",
        names=["w", "d", "sigma", "omega", "type", "run"],
    )

    # Create plot
    plot_freq(df)


if __name__ == "__main__":
    main()
