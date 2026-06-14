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
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, WhiteKernel
from sklearn.gaussian_process.kernels import Matern
from sklearn.preprocessing import StandardScaler


def computeAmplitude(dataFile: str, doLoo: bool) -> Tuple[np.ndarray, np.ndarray]:
    d, w = extract_width_depth(dataFile)
    print(colors.OKBLUE + f"Processing d = {d:.2f}, w = {w:.2f}" + colors.ENDC)
    try:
        x, y, data = readFieldsBySection(dataFile)
    except Exception as e:
        print(colors.FAIL + f"Error reading data from {dataFile}: {e}" + colors.ENDC)
        return np.array([]), np.array([])

    rms = getRMS(data)

    if doLoo:
        Loo = np.max(np.abs(rms), axis=1)
        return x, Loo
    else:
        # integrate from 0 to ymax
        ymax = 150
        yindx = np.argmin(np.abs(y - ymax))

        L2 = np.sqrt(np.trapezoid(rms[:, :yindx] ** 2, y[:yindx]))

        return x, L2

def computeNx(dataFile: str, doLoo: bool) -> Tuple[np.ndarray, np.ndarray]:
    x, A = computeAmplitude(dataFile, doLoo)
    if len(x) == 0 or len(A) == 0:
        return np.array([]), np.array([])
    

    # filter indices of x such that x <=0
    idx = np.where(x <= 0)[0]

    A0_arg = np.argmin(A[idx])
    A0 = A[A0_arg]
    A = A[A0_arg:]
    x = x[A0_arg:]

    Nx = np.log(A / A0)

    return x, Nx 

def computeDeltaN(
    w: float, x: np.ndarray, Nx: np.ndarray, x_flat: np.ndarray, Nx_flat: np.ndarray, x_start: float = 50, x_end: float = 250
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


class GPRResult:
    """Container for GPR interpolation results."""
    def __init__(
        self,
        X_grid: np.ndarray,
        Y_grid: np.ndarray,
        Z_grid: np.ndarray,
        sigma_grid: np.ndarray,
        widths: np.ndarray,
        depths: np.ndarray,
        deltaNx: np.ndarray,
        gp: GaussianProcessRegressor,
        scaler: StandardScaler = None,
    ):
        self.X_grid = X_grid
        self.Y_grid = Y_grid
        self.Z_grid = Z_grid
        self.sigma_grid = sigma_grid
        self.widths = widths
        self.depths = depths
        self.deltaNx = deltaNx
        self.gp = gp
        self.scaler = scaler


def gpr_interpolate(
    depths: np.ndarray,
    widths: np.ndarray,
    deltaNx: np.ndarray,
    w_max: float = 90.0,
    grid_size: int = 200,
) -> GPRResult:
    """
    Fits a GPR model and interpolates/extrapolates deltaN values over a grid.
    
    Returns a GPRResult containing grids and the fitted model.
    """
    # Make copies to avoid modifying input arrays
    depths = depths.copy()
    widths = widths.copy()
    deltaNx = deltaNx.copy()

    # Remove NaN/invalid values
    valid_mask = ~np.isnan(deltaNx) & ~np.isinf(deltaNx)
    depths = depths[valid_mask]
    widths = widths[valid_mask]
    deltaNx = deltaNx[valid_mask]

    # Store original data points for plotting
    orig_widths = widths.copy()
    orig_depths = depths.copy()
    orig_deltaNx = deltaNx.copy()

    # Add denser boundary conditions (zero values at edges)
    n_boundary = 15
    max_depth = np.max(depths)
    max_width = w_max

    # Along w=0 boundary
    for d in np.linspace(0, max_depth, n_boundary):
        depths = np.append(depths, d)
        widths = np.append(widths, 0)
        deltaNx = np.append(deltaNx, 0)

    # Along d=0 boundary
    for w in np.linspace(0, max_width, n_boundary):
        depths = np.append(depths, 0)
        widths = np.append(widths, w)
        deltaNx = np.append(deltaNx, 0)

    # Corner point (0,0)
    depths = np.append(depths, 0)
    widths = np.append(widths, 0)
    deltaNx = np.append(deltaNx, 0)

    # Add soft prior points from empirical model in sparse regions
    # deltaN_approx = 0.1 * w * tanh(44 * d / w) with added noise
    prior_w = np.array([60, 70, 80, 50, 60, 70, 80])
    prior_d = np.array([2.0, 2.0, 2.0, 2.5, 2.5, 2.5, 2.5])
    prior_deltaN = 0.1 * prior_w * np.tanh(44 * prior_d / prior_w)
    
    # Add these soft priors
    for w, d, dN in zip(prior_w, prior_d, prior_deltaN):
        widths = np.append(widths, w)
        depths = np.append(depths, d)
        deltaNx = np.append(deltaNx, dN)

    # Training data
    X_train = np.column_stack((widths, depths))
    y_train = deltaNx

    # Feature scaling for better GPR performance
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)

    # Kernel: Matérn kernel (better for physical processes) with anisotropic lengthscales
    kernel = (
        Matern(
            length_scale=[1.0, 1.0],  # After scaling, start with unit lengthscales
            length_scale_bounds=[(1e-2, 1e2), (1e-2, 1e2)],
            nu=2.5
        ) +
        WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-10, 1e-1))
    )

    gp = GaussianProcessRegressor(
        kernel=kernel,
        n_restarts_optimizer=20,
        normalize_y=True,
        random_state=42
    )

    gp.fit(X_train_scaled, y_train)

    print(colors.OKBLUE + f"Learned kernel: {gp.kernel_}" + colors.ENDC)
    print(colors.OKBLUE + f"Log-marginal-likelihood: {gp.log_marginal_likelihood_value_:.3f}" + colors.ENDC)

    # Prediction grid
    X = np.linspace(0, w_max, grid_size)
    Y = np.linspace(0, np.max(orig_depths), grid_size)
    X_grid, Y_grid = np.meshgrid(X, Y)
    XY = np.column_stack((X_grid.ravel(), Y_grid.ravel()))

    # Scale prediction points using the same scaler
    XY_scaled = scaler.transform(XY)

    Z_pred, sigma = gp.predict(XY_scaled, return_std=True)
    Z_grid = Z_pred.reshape(X_grid.shape)
    sigma_grid = sigma.reshape(X_grid.shape)

    # Clip negative predictions (deltaN should be non-negative)
    Z_grid = np.clip(Z_grid, 0, None)

    return GPRResult(
        X_grid, Y_grid, Z_grid, sigma_grid,
        orig_widths, orig_depths, orig_deltaNx,
        gp, scaler
    )


def write_deltaN_file(result: GPRResult, filename: str) -> None:
    """Writes the interpolated deltaN grid to a .dat file."""
    X = result.X_grid[0, :]
    Y = result.Y_grid[:, 0]
    Z_grid = result.Z_grid

    with open(filename, "w") as f:
        for i in range(len(Y)):
            if i == 0:
                f.write(f"{0:.6f} {0:.6f} {0:.6f}\n")
                for j in range(len(X)):
                    f.write(f"{X[j]:.6f} {0:.6f} {0:.6f}\n")
                f.write("\n")

            for j in range(len(X)):
                if j == 0:
                    f.write(f"{0:.6f} {Y[i]:.6f} {0:.6f}\n")
                f.write(f"{X[j]:.6f} {Y[i]:.6f} {Z_grid[i, j]:.6f}\n")
            f.write("\n")

    print(colors.OKGREEN + f"Interpolated data written to {filename}" + colors.ENDC)


def plot_deltaN_grid(result: GPRResult, ax: plt.Axes = None) -> plt.Figure:
    """Plots the interpolated deltaN grid with scatter points and contours."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.get_figure()

    c = ax.pcolormesh(result.X_grid, result.Y_grid, result.Z_grid, shading="auto", cmap="viridis")
    ax.scatter(result.widths, result.depths, c=result.deltaNx, edgecolors="k", cmap="viridis", s=100)
    
    for w, d, val in zip(result.widths, result.depths, result.deltaNx):
        ax.text(w, d - 0.05, f"{val:.2f}", color="white", ha="center", va="bottom", fontsize=8)

    contour_levels = [1, 2, 3, 4, 5]
    CS = ax.contour(result.X_grid, result.Y_grid, result.Z_grid, levels=contour_levels, colors="white", linewidths=1.2)
    ax.clabel(CS, inline=True, fontsize=8, fmt="%d")

    fig.colorbar(c, ax=ax, label="Delta N")
    ax.set_xlabel("Width (w)")
    ax.set_ylabel("Depth (d)")
    ax.set_title("Interpolated and Extrapolated Delta N")

    return fig




