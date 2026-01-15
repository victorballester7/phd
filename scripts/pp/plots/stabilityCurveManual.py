import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from pp.colors import colors
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

plt.style.use("plots/style/customvictor.mplstyle")


class runType:
    def __init__(self, codeDatFile, label, weightX, endFrac=1.0):
        self.codeDatFile = codeDatFile
        self.label = label
        self.weightX = weightX
        self.endFrac = endFrac


inc2d = runType("inc2d", "Inc. 2d", 0.0333, 1.0)
ma02_2d = runType("Ma0.2_2d", "Ma = 0.2 2d", 0.0333)
ma04_2d = runType("Ma0.4_2d", "Ma = 0.4 2d", 0.0333)
ma06_2d = runType("Ma0.6_2d", "Ma = 0.6 2d", 0.035, 0.95)
ma08_2d = runType("Ma0.8_2d", "Ma = 0.8 2d", 0.028, 0.48)
inc3d = runType("inc3d", "Inc. 3d", 0.05, 0.58)

PLOT = [inc2d, ma02_2d, ma04_2d, ma06_2d, ma08_2d]

COLORS = plt.cm.tab10.colors  # or any colormap you like
alpha = 0.15

PLOT_SCATTER = True  # Whether to plot scatter points


def plot_stability_diagram(df, output_path):
    """
    Plot stability diagram and extract bifurcation line using SVM.
    """
    fig, ax = plt.subplots(figsize=(7, 5))

    legend_handles = []

    for i, run in enumerate(PLOT):
        # --------------------------------------------------
        # 1. Prepare data
        # --------------------------------------------------

        # Map type -> binary labels

        dftmp = df.loc[df["run"] == run.codeDatFile].copy()

        # Limit to d>= 2.5
        # dftmp = dftmp.loc[dftmp["d"] >= 2.25]

        dftmp["label"] = (dftmp["type"] != "equilibrium").astype(int)

        X = dftmp[["w", "d"]].values
        y = dftmp["label"].values

        # --------------------------------------------------
        # 2. Build classifier pipeline
        # --------------------------------------------------

        # Standardization is essential due to different scales
        model = Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "svm",
                    SVC(
                        kernel="rbf",  # smooth nonlinear boundary
                        C=500.0,  # softness of margin
                        gamma="scale",  # good default
                    ),
                ),
            ]
        )

        model.fit(X, y)

        # --------------------------------------------------
        # 3. Create grid for decision function
        # --------------------------------------------------

        w_min, w_max = X[:, 0].min(), X[:, 0].max()
        d_min, d_max = X[:, 1].min(), X[:, 1].max()

        w_grid = np.linspace(w_min, w_max, 400)
        d_grid = np.linspace(d_min, d_max, 400)
        W, D = np.meshgrid(w_grid, d_grid)

        grid_points = np.c_[W.ravel(), D.ravel()]

        decision_values = model.decision_function(grid_points)
        decision_values = decision_values.reshape(W.shape)

        # --------------------------------------------------
        # 4. Plot
        # --------------------------------------------------

        # Stable points
        stable = dftmp[dftmp["label"] == 0]
        unstable = dftmp[dftmp["label"] == 1]

        color = COLORS[i % len(COLORS)]

        if PLOT_SCATTER:
            ax.scatter(
                stable["w"],
                stable["d"],
                marker="o",
                label="Stable",
                color=color,
                alpha=alpha,
            )
            ax.scatter(
                unstable["w"],
                unstable["d"],
                marker="x",
                label="Unstable",
                color=color,
                alpha=alpha,
            )

        # Bifurcation line: decision = 0
        ax.contour(
            W,
            D,
            decision_values,
            levels=[0],
            colors=color,
            linewidths=2,
            linestyles="-",
        )

        # ---- manual legend entry for this bifurcation line
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                lw=2,
                label=run.label,  # or run.name / run.codeDatFile
            )
        )

    ax.set_xlabel(r"$w/\delta^*$")
    ax.set_ylabel(r"$d/\delta^*$")
    # ax.set_xlim(0, df["w"].max() * 1.005)
    # ax.set_ylim(0, df["d"].max() * 1.05)
    legend_handles = [
        Line2D([0], [0], marker="o", color="k", linestyle="None", label="Equilibrium"),
        Line2D([0], [0], marker="x", color="k", linestyle="None", label="Unstable"),
    ] + legend_handles

    ax.legend(
        handles=legend_handles,
        loc=(1.02, 0.0),
        frameon=True,
    )

    ax.grid(True)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def main():
    """Main execution function."""
    # Set up paths
    script_path = os.path.dirname(os.path.abspath(__file__))
    data_filename = os.path.join(script_path, "../../../data/stability/stability.dat")
    output_path = os.path.join(script_path, "../../../images/stabilityCurveNew.pdf")

    # Load data
    data_frame = pd.read_csv(
        data_filename,
        sep=" ",
        comment="#",
        names=["w", "d", "sigma", "omega", "type", "run"],
    )

    # Create plot
    plot_stability_diagram(data_frame, output_path)


if __name__ == "__main__":
    main()
