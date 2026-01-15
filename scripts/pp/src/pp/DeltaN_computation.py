import numpy as np
from pp.colors import colors
from scipy.interpolate import (
    griddata,
    NearestNDInterpolator,
    CloughTocher2DInterpolator,
)
from scipy.ndimage import gaussian_filter, median_filter
import matplotlib.pyplot as plt
from pp.fileManagement import extract_width_depth, readFieldsBySection
from pp.filterData import getRMS
from pp.colors import colors
from typing import Tuple
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C



def computeNx(dataFile: str, doLoo: bool) -> Tuple[float, float, np.ndarray, np.ndarray]:
    d, w = extract_width_depth(dataFile)
    print(colors.OKBLUE + f"Processing d = {d:.2f}, w = {w:.2f}" + colors.ENDC)
    x, y, data = readFieldsBySection(dataFile)

    rms = getRMS(data)

    if doLoo:
        Loo = np.max(np.abs(rms), axis=1)
        A = Loo
    else:
        # integrate from 0 to ymax
        ymax = 150
        yindx = np.argmin(np.abs(y - ymax))

        L2 = np.sqrt(np.trapezoid(rms[:, :yindx] ** 2, y[:yindx]))

        A = L2

    A0_arg = np.argmin(A)
    A0 = A[A0_arg]
    A = A[A0_arg:]
    x = x[A0_arg:]

    Nx = np.log(A / A0)

    return d, w, x, Nx 

def computeDeltaN(
    w: float, x: np.ndarray, Nx: np.ndarray, x_flat: np.ndarray, Nx_flat: np.ndarray
) -> float:
    """
    Computes the difference in N factor between a given case and a flat plate case at in a prescribed window of x values.
    """

    # interpolate x_flat, Nx_flat based on the values of x
    Nx_flat_interp = np.interp(x, x_flat, Nx_flat)

    # compute the average of deltaNx in the window w+50 < x < w+250
    x_start = w + 50
    x_end = w + 250
    indices = np.where((x >= x_start) & (x <= x_end))[0]
    if len(indices) == 0:
        print(
            colors.WARNING
                + f"No points found in the window {x_start} < x < {x_end} for w = {w}"
                + colors.ENDC
        )
        return np.nan

    deltaNx = Nx - Nx_flat_interp
    deltaNx_avg = np.mean(deltaNx[indices])
    print(colors.OKGREEN + f"Computed deltaN for w = {w}: {deltaNx_avg}" + colors.ENDC)

    return deltaNx_avg.astype(float)


def interpolate_extrapolate(
    depths: np.ndarray, widths: np.ndarray, deltaNx: np.ndarray, filename: str
) -> None:
    """
    Interpolates and extrapolates the deltaN values over a grid defined by the maximum and minimum values of d and w.
    """

    # add 0 values of depth and width
    for d in np.unique(depths):
        depths = np.append(depths, d)
        widths = np.append(widths, 0)
        deltaNx = np.append(deltaNx, 0)
    for w in np.unique(widths):
        depths = np.append(depths, 0)
        widths = np.append(widths, w)
        deltaNx = np.append(deltaNx, 0)
    depths = np.append(depths, 0)
    widths = np.append(widths, 0)
    deltaNx = np.append(deltaNx, 0)


    # training data
    X_train = np.column_stack((widths, depths))
    y_train = deltaNx

    # kernel: constant * RBF
    kernel = RBF([1, 10], [(1e-2, 1e2),(1e-1,1e3)])  # anisotropic lengthscales
    # kernel = C(1.0, (1e-3, 1e3)) * RBF([10, 10], (1e-2, 1e2))  # anisotropic lengthscales
    gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10, alpha=1e-6, normalize_y=True)

    gp.fit(X_train, y_train)

    # prediction grid
    sizeGrid = 100
    X = np.linspace(np.min(widths), 70, sizeGrid)
    Y = np.linspace(np.min(depths), np.max(depths), sizeGrid)
    X_grid, Y_grid = np.meshgrid(X, Y)
    XY = np.column_stack((X_grid.ravel(), Y_grid.ravel()))

    Z_pred, sigma = gp.predict(XY, return_std=True)
    Z_grid = Z_pred.reshape(X_grid.shape)

    # sizeGrid = 100
    # X = np.linspace(np.min(widths), 70, sizeGrid)
    # # X = np.linspace(np.min(widths), np.max(widths), sizeGrid)
    # Y = np.linspace(np.min(depths), np.max(depths), sizeGrid)
    # X_grid, Y_grid = np.meshgrid(X, Y)

    # # Interpolate/extrapolate over the grid
    # points = np.column_stack((widths, depths))

    # Z_grid = griddata(
    #     points, deltaNx, (X_grid, Y_grid), method="linear", fill_value=np.nan
    # )

    # # Step 2: Fill NaNs using nearest neighbor extrapolation
    # nearest = NearestNDInterpolator(points, deltaNx)
    # Z_grid = np.where(np.isnan(Z_grid), nearest(X_grid, Y_grid), Z_grid)

    # Z_grid = gaussian_filter(Z_grid, sigma=2)

    # Write to .dat file
    with open(filename, "w") as f:
        for i in range(len(Y)):
            if i == 0:
                # Write the border points
                f.write(f"{0:.6f} {0:.6f} {0:.6f}\n")
                for j in range(len(X)):
                    f.write(f"{X[j]:.6f} {0:.6f} {0:.6f}\n")
                f.write("\n")  # Newline to separate rows

            for j in range(len(X)):
                if j == 0:
                    # Write the border points
                    f.write(f"{0:.6f} {Y[i]:.6f} {0:.6f}\n")
                x = X[j]
                y = Y[i]
                z = Z_grid[i, j]
                # if np.isnan(z):
                #     continue  # Skip undefined regions
                f.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
            f.write("\n")  # Newline to separate rows

    print(
        colors.OKGREEN
        + f"Interpolated and extrapolated data written to {filename}"
        + colors.ENDC
    )

    fig, ax = plt.subplots(figsize=(8, 6))
    c = ax.pcolormesh(X_grid, Y_grid, Z_grid, shading="auto", cmap="viridis")
    # plot a ball with the value as label text above the vall of deltaN at each (w,d) point
    sc = ax.scatter(widths, depths, c=deltaNx, edgecolors="k", cmap="viridis", s=100)
    for (i, j, val) in zip(widths, depths, deltaNx):
        ax.text(i, j - 0.05, f"{val:.2f}", color="white", ha="center", va="bottom", fontsize=8)

    # add contour lines for Delta N = 1, 2, 3, 4, 5
    contour_levels = [1, 2, 3, 4, 5]
    CS = ax.contour(X_grid, Y_grid, Z_grid, levels=contour_levels, colors="white", linewidths=1.2)
    ax.clabel(CS, inline=True, fontsize=8, fmt="%d")  # label contours

    fig.colorbar(c, ax=ax, label="Delta N")
    ax.set_xlabel("Width (w)")
    ax.set_ylabel("Depth (d)")
    ax.set_title("Interpolated and Extrapolated Delta N")
