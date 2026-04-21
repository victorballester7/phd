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
from scipy.optimize import minimize


plt.style.use("plots/style/customvictor.mplstyle")


class runType:
    def __init__(self, codeDatFile, label, weightX, endFrac=1.0):
        self.codeDatFile = codeDatFile
        self.label = label
        self.weightX = weightX
        self.endFrac = endFrac


incRe10002d = runType("inc2dRe1000", "Inc. Re = 1000 2d", 0.0333, 1.0)
incRe20002d = runType("inc2dRe2000", "Inc. Re = 2000 2d", 0.0333, 1.0)
ma005Re1000_2d = runType("Ma0.05Re1000_2d", "Ma = 0.05 Re = 1000 2d", 1.0)
ma02Re1000_2d = runType("Ma0.2Re1000_2d", "Ma = 0.2 Re = 1000 2d", 0.0333)
ma04Re1000_2d = runType("Ma0.4Re1000_2d", "Ma = 0.4 Re = 1000 2d", 0.0333)
ma06Re1000_2d = runType("Ma0.6Re1000_2d", "Ma = 0.6 Re = 1000 2d", 0.035, 0.95)
ma08Re1000_2d = runType("Ma0.8Re1000_2d", "Ma = 0.8 Re = 1000 2d", 0.028, 0.48)
incRe10003d = runType("inc3dRe1000", "Inc. Re = 1000 3d", 0.05, 0.58)

PLOT = [
    incRe10002d,
    # incRe20002d,
    # incRe10003d,
    ma005Re1000_2d,
    ma02Re1000_2d,
    ma04Re1000_2d,
    ma06Re1000_2d,
    ma08Re1000_2d,
]


COLORS = plt.cm.tab10.colors  # or any colormap you like
alpha = 0.15

PLOT_SCATTER = True  # Whether to plot scatter points


def scaled_norm(v, l_w, l_d):
    return np.sqrt((v[0] / l_w) ** 2 + (v[1] / l_d) ** 2)


def scaled_dist(p, q, l_w, l_d):
    return scaled_norm(p - q, l_w, l_d)


def dist_to_set(point, cloud, l_w, l_d):
    diffs = cloud - point
    dists = np.sqrt((diffs[:, 0] / l_w) ** 2 + (diffs[:, 1] / l_d) ** 2)
    return np.min(dists)


# def dist_to_set(point, cloud, tau=0.05):
#     diffs = cloud - point
#     dists = np.sqrt((diffs[:, 0] / L_W) ** 2 + (diffs[:, 1] / L_D) ** 2)
#     return -tau * np.log(np.sum(np.exp(-dists / tau)) + 1e-8)  # add small term to avoid log(0)


def energy(X, A, B, N, l_w, l_d, lam=0.2, mu=0.001, nu=0.1):
    X = X.reshape(N, 2)

    E = 0.0
    tmp = []

    # data fidelity
    alpha = 0.5
    weight_distance = 2.0
    weight_smooth = 0.1
    bigPenalty = 1
    for i in range(N):
        dA = dist_to_set(X[i], A, l_w, l_d)
        dB = dist_to_set(X[i], B, l_w, l_d)
        # weight = 1.0 / (1.0 + 0.01 * X[i,0])
        E += weight_distance * (alpha * dA**2 + (1 - alpha) * dB**2)

        # impose positive gradient in the w direction
        # if i > 0:
        #     dw = X[i, 0] - X[i - 1, 0]
        #     if dw < 0:
        #         E += penalty * dw ** 2

    tmp.append(E)
    # concAB = np.concatenate((A, B), axis=0)
    # max_d = np.max(concAB[:, 1])
    # max_w = np.max(concAB[:, 0])
    # print("max d, w:", max_d, max_w)

    # impose start point for d and end point for w
    # E += penalty * (X[0, 1] - max_d) ** 2
    # E += penalty * (X[-1, 0] - max_w) ** 2

    # length regularization
    seg_lengths = np.array(
        [scaled_dist(X[i], X[i + 1], l_w, l_d) for i in range(N - 1)]
    )

    # compute for each point X[i] the mean distance to all the other points in the curve.
    # total_distances =np.zeros((N, N))
    # for i in range(N):
    #     for j in range(i + 1, N):
    #         d = scaled_dist(X[i], X[j])
    #         total_distances[i, j] = d
    #         total_distances[j, i] = d
    # mean_distances = np.mean(total_distances, axis=1)
    # E += 0.4 * np.sum(mean_distances**2)
    # print("Segment lengths:", seg_lengths)

    # segments smaller than eps are penalized quadratically
    # eps = 0.5
    # E += bigPenalty * np.sum(np.minimum(0, seg_lengths - eps) ** 2)

    # impose first point is a t d=4
    if np.abs(X[0, 1] - 4) > 0.1: 
        E += bigPenalty * (X[0, 1] - 4) ** 2

    # curve_min_dist_to_itself = np.array([
    #     np.min([scaled_dist(X[i], X[j], l_w, l_d) for j in range(N) if j != i]) for i in range(N)
    # ])

    # # impose minimum distance between points to avoid self-intersection
    # E+= bigPenalty * np.sum(np.minimum(0, curve_min_dist_to_itself - 0.3) ** 2)

    # gradient penalization
    dX = np.gradient(X, axis=0)
    # encourage positive gradient in w direction
    dX_w = dX[:, 0]
    E += 0.6 * np.sum(np.minimum(0, dX_w) ** 2)

    tmp.append(E - tmp[-1])

    # # equal spacing penalty
    # # ell = np.mean(seg_lengths)
    # # E += 0.4 * np.sum((seg_lengths - ell) ** 2)

    # smoothness regularization
    for i in range(1, N - 1):
        E += weight_smooth * scaled_norm(X[i + 1] - 2 * X[i] + X[i - 1], l_w, l_d) ** 2
    tmp.append(E - tmp[-1])

    # curvature regularization
    # ddX = np.gradient(dX, axis=0)
    # ddX_norms = np.array([scaled_norm(ddX[i],l_w,l_d) for i in range(N)])
    # E += 0.4 * np.sum(ddX_norms**2)
    tmp.append(E - tmp[-1])

    # print(f"Energy components: data={tmp[0]:.4f}, grad={tmp[1]:.4f}, smooth={tmp[2]:.4f}, curvature={tmp[3]:.4f}")

    return E


def getBifCurve(A, B, N, l_w, l_d):
    """
    Find boundary curve between equilibrium (A) and non-equilibrium (B) points.

    Strategy:
    1. Initialize along a simple boundary
    2. Minimize energy with soft constraints
    3. Use moderate N (10-20 points)
    """

    # ===== SMART INITIALIZATION =====
    # Find boundary points: max depth for equilibrium, min depth for non-equilibrium
    w_grid = np.linspace(A[:, 0].min(), A[:, 0].max(), N)
    d_init = np.zeros(N)

    for i, w in enumerate(w_grid):
        # For this w, find nearby equilibrium and non-equilibrium points
        nearby_A = A[np.abs(A[:, 0] - w) < 20]  # within width 20
        nearby_B = B[np.abs(B[:, 0] - w) < 20]

        if len(nearby_A) > 0 and len(nearby_B) > 0:
            # Initialize midway between max equilibrium depth and min non-eq depth
            d_init[i] = 0.5 * (nearby_A[:, 1].max() + nearby_B[:, 1].min())
        elif len(nearby_A) > 0:
            d_init[i] = nearby_A[:, 1].max()
        elif len(nearby_B) > 0:
            d_init[i] = nearby_B[:, 1].min()
        else:
            # Linear interpolation fallback
            d_init[i] = 4.0 - (w - w_grid[0]) * 3.0 / (w_grid[-1] - w_grid[0])

    X0 = np.column_stack([w_grid, d_init]).flatten()

    # print(f"Initial energy: {energy_simple(X0, A, B, N):.2f}")
    print(f"Initial energy: {energy(X0, A, B, N, l_w, l_d):.2f}")

    # ===== OPTIMIZATION =====
    # Set bounds to keep points in reasonable range
    w_min, w_max = A[:, 0].min() - 5, max(A[:, 0].max(), B[:, 0].max()) + 5
    d_min, d_max = 0.5, max(A[:, 1].max(), B[:, 1].max()) + 0.5

    bounds = [(w_min, w_max), (d_min, d_max)] * N
    result = minimize(
        energy,
        X0,
        args=(A, B, N, l_w, l_d),
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 2000},
    )

    print(f"Optimization success: {result.success}")
    print(f"Final energy: {result.fun:.2f}")

    gamma = result.x.reshape(N, 2)

    # ===== POST-PROCESS: Sort by w to ensure monotonicity =====
    # sorted_idx = np.argsort(gamma[:, 0])
    # gamma = gamma[sorted_idx]

    return gamma


def plot_stability_diagram(df, output_path):
    """
    Plot stability diagram and extract bifurcation line using SVM.
    """
    # 3 x 2 figures
    rows, cols = 3, int(np.ceil(len(PLOT) / 3))
    _, ax = plt.subplots(
        rows, cols, figsize=(5 * cols, 4.3 * rows), sharex=True, sharey=True
    )

    ax = ax.flatten()

    curves = []

    for i, run in enumerate(PLOT):
        # --------------------------------------------------
        # 1. Prepare data
        # --------------------------------------------------

        # Map type -> binary labels

        dftmp = df.loc[df["run"] == run.codeDatFile].copy()

        # Limit to d>= 2.5
        # dftmp = dftmp.loc[dftmp["d"] >= 2.5]

        # dftmp["label"] = (dftmp["type"] != "equilibrium").astype(int)

        # X = dftmp[["w", "d"]].values
        # y = dftmp["label"].values

        A = dftmp.loc[
            (dftmp["type"] == "equilibrium") | (dftmp["type"] == "bursts"),
            ["w", "d"],
        ].values

        B = dftmp.loc[
            (dftmp["type"] != "equilibrium") & (dftmp["type"] != "bursts"),
            ["w", "d"],
        ].values

        concatAB = np.concatenate((A, B), axis=0)

        std_w = np.std(concatAB[:, 0])
        std_d = np.std(concatAB[:, 1])
        print(f"std w: {std_w:.2f}, std d: {std_d:.2f}")
        print(f"Running {run.label}...")
        l_w = 80
        l_d = 4
        N = 20
        gamma = getBifCurve(A, B, N, l_w, l_d)
        # plot as np array to be easily copied
        print("gamma points")
        print("np.array([")
        for point in gamma:
            print(f"    [{point[0]:.2f}, {point[1]:.2f}],")
        print("])")

        color = COLORS[i % len(COLORS)]

        # gradient in both directions of gamma

        # tangents = np.gradient(gamma, axis=0)
        # print(tangents)

        if PLOT_SCATTER:
            ax[i].scatter(
                A[:, 0],
                A[:, 1],
                marker="o",
                label="Stable",
                color=color,
                alpha=alpha,
            )
            ax[i].scatter(
                B[:, 0],
                B[:, 1],
                marker="x",
                label="Unstable",
                color=color,
                alpha=alpha,
            )

        ax[i].plot(
            gamma[:, 0],
            gamma[:, 1],
            "*-",
            color=color,
            linewidth=0.5,
            markersize=4,
            label=run.label,
        )

        curves.append([gamma, run.label, color])

        ax[i].set_title(run.label)
        ax[i].set_xlabel(r"$w/\delta^*$")
        ax[i].set_ylabel(r"$d/\delta^*$")

        ax[i].grid(True)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

    _, ax = plt.subplots(figsize=(6, 5))
    for gamma, label, color in curves:
        ax.plot(gamma[:, 0], gamma[:, 1], color=color, linewidth=2, label=label)

    # tmp #################################
    D_points = np.array([4, 4, 2, 2, 1.5])
    W_points = np.array([14, 19, 24, 40 , 80])

    ax.scatter(W_points, D_points, color="black", marker="x", label="Full 3d simulations")
    #######################################

    ax.set_xlabel(r"$w/\delta^*$")
    ax.set_ylabel(r"$d/\delta^*$")
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.savefig(
        os.path.join(os.path.dirname(output_path), "stabilityCurvesCombined.pdf")
    )
    plt.close()


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
