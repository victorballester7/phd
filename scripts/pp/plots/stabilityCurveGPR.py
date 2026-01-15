import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from pp.colors import colors
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn.preprocessing import StandardScaler

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

PLOT = [inc2d]

COLORS = plt.cm.tab10.colors  # or any colormap you like

PLOT_SCATTER = True  # Whether to plot scatter points

def dist(set1, set2, weight_x, weight_y):
    """
    For each point in set1, compute the squared distance
    to the nearest point in set2 (anisotropic metric).
    """

    dx = weight_x * (set1[:, None, 0] - set2[None, :, 0])
    dy = weight_y * (set1[:, None, 1] - set2[None, :, 1])


    dist = np.min(dx**2 + dy**2, axis=1)  # shape (N1,)


    # create softmin distance
    
    # dist = - (1/beta) * np.log(np.sum(np.exp(-beta * (dx**2 + dy**2)), axis=1))


    return np.hstack((set1, dist[:, None]))


def gpr(df, run_label, stable_type, weight_x=0.033333, weight_y=1.0):
    stable = df.loc[
        (df["run"] == run_label) & (df["type"] == stable_type), ["w", "d"]
    ].to_numpy()
    unstable = df.loc[
        (df["run"] == run_label) & (df["type"] != stable_type), ["w", "d"]
    ].to_numpy()

    distsS2US = dist(stable, unstable, weight_x, weight_y)
    distsS2US[:, 2] *= -1  # negative distances for stable points
    distsUS2S = dist(unstable, stable, weight_x, weight_y)

    data = np.vstack((distsS2US, distsUS2S))

    X = data[:, :2]
    y = data[:, 2]
    scaler = StandardScaler()
    Xn = scaler.fit_transform(X)

    kernel = 1.0 * RBF(
        length_scale=[1.0, 1.0], length_scale_bounds=(1e-2, 1e2)
    ) + WhiteKernel(noise_level=1e-6)

    gpr = GaussianProcessRegressor(
        kernel=kernel,
        normalize_y=True,
        n_restarts_optimizer=5,
    )

    gpr.fit(Xn, y)

    # grid in original coordinates
    w = np.linspace(df["w"].min(), df["w"].max(), 300)
    d = np.linspace(df["d"].min(), df["d"].max(), 300)
    W, D = np.meshgrid(w, d)

    Xg = np.column_stack([W.ravel(), D.ravel()])
    Xg_n = scaler.transform(Xg)

    yg, yg_std = gpr.predict(Xg_n, return_std=True)
    Z = yg.reshape(W.shape)

    return W, D, Z, stable, unstable


def plot_stability_diagram(df, save_path):
    fig, ax = plt.subplots(figsize=(8, 6))

    for i, run in enumerate(PLOT):
        run_label = run.codeDatFile
        print(colors.OKBLUE + f"Processing run: {run_label}" + colors.ENDC)

        W, D, Z, stable, unstable = gpr(df, run_label, "equilibrium", run.weightX, 1.0)

        color = COLORS[i % len(COLORS)]

        if PLOT_SCATTER:
            ax.scatter(
                stable[:, 0],
                stable[:, 1],
                marker="o",
                label="Stable",
                color=color,
                alpha=0.15,
            )
            ax.scatter(
                unstable[:, 0],
                unstable[:, 1],
                marker="x",
                label="Unstable",
                color=color,
                alpha=0.15,
            )

        contour_set = ax.contour(
            W,
            D,
            Z,
            levels=[0.0],
            alpha=0.0,
        )

        # Extract and plot partial contour paths
        for level_paths in contour_set.allsegs:
            # print(level_paths)
            for path_vertices in level_paths:
                # print(path_vertices)
                n_points = len(path_vertices)

                # Calculate indices for the desired fraction
                start_idx = 0
                end_idx = int(run.endFrac * n_points)

                # Extract and plot the partial path
                partial_vertices = path_vertices[start_idx:end_idx]
                ax.plot(
                    partial_vertices[:, 0],
                    partial_vertices[:, 1],
                    linewidth=2,
                    color=color,
                    label=run.label,
                )

    ax.set_xlabel(r"$w/\delta^*$")
    ax.set_ylabel(r"$d/\delta^*$", rotation=0, labelpad=20)
    ax.set_xlim(0, df["w"].max() * 1.005)
    ax.set_ylim(0, df["d"].max() * 1.05)

    ax.legend()

    plt.tight_layout()
    plt.savefig(save_path)
    print(colors.OKGREEN + f"✓ Plot saved to: {save_path}" + colors.ENDC)


def main():
    """Main execution function."""
    # Set up paths
    script_path = os.path.dirname(os.path.abspath(__file__))
    data_filename = os.path.join(script_path, "../../../data/stability/stability.dat")
    output_path = os.path.join(script_path, "../../../images/stabilityCurve.pdf")

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
