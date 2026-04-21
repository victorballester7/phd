import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

MARKERS = ["o", "s", "^", "P", "*", "X", "h", "8"]

runs = ["inc2dRe1000"]

depth_same_color = False


def empiricalFormula(wOverd, M, n) -> float:
    return n / (1 + M + 1.0 / np.tanh(np.pi / wOverd))


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

        widthOverDepth = []
        fLOverU = []

        M = 0
        n = 1

        # label_done = False
        for i, depth in enumerate(depths):
            df_depth = df_run.loc[df_run["d"] == depth]

            for w in df_depth["w"]:
                widthOverDepth.append(w / depth)
                omega = df_depth.loc[df_depth["w"] == w, "omega"].values[0]
                st = omega / 2 / np.pi * w
                fLOverU.append(st)

        widthOverDepth = np.array(widthOverDepth)

        withOverDepth_linspace = np.linspace(
            np.min(widthOverDepth), np.max(widthOverDepth), 100
        )
        fLOverU_empirical = empiricalFormula(withOverDepth_linspace, M, n)

        ax.plot( widthOverDepth, fLOverU, label=f"Run: {run}", marker='o')
        ax.plot( withOverDepth_linspace, fLOverU_empirical, label=f"Empirical Formula", marker='x')



    ax.set_xlabel("w/d")
    ax.set_ylabel("St = fL/U")
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
