import numpy as np
from pp.colors import colors
from scipy.ndimage import gaussian_filter, median_filter
import matplotlib.pyplot as plt
from pp.fileManagement import extract_depth_width, readFieldsBySection
from pp.filterData import getRMS
from typing import Tuple
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, WhiteKernel
from sklearn.gaussian_process.kernels import Matern, Kernel, Hyperparameter
from sklearn.preprocessing import StandardScaler


def computeAmplitude(dataFile: str, doLoo: bool, field: str) -> Tuple[np.ndarray, np.ndarray]:
    d, w = extract_depth_width(dataFile)
    print(colors.OKBLUE + f"Processing d = {d:.2f}, w = {w:.2f}" + colors.ENDC)
    try:
        x, y, data = readFieldsBySection(dataFile)
    except Exception as e:
        print(colors.FAIL + f"Error reading data from {dataFile}: {e}" + colors.ENDC)
        return np.array([]), np.array([])

    f = np.zeros_like(data[:, :, 0])  # Initialize f with the same shape as one field component
    
    match field:
        case "u":
            f = data[:, :, 2]
        case "v":
            f = data[:, :, 3]
        case "|u|":
            f = np.sqrt(data[:, :, 2] ** 2 + data[:, :, 3] ** 2)
        case "rms":
            f = getRMS(data)


    if doLoo:
        Loo = np.max(np.abs(f), axis=1)
        return x, Loo
    else:
        # integrate from 0 to ymax
        ymax = 150
        yindx = np.argmin(np.abs(y - ymax))

        L2 = np.sqrt(np.trapezoid(f[:, :yindx] ** 2, y[:yindx]))

        return x, L2

def computeNx(dataFile: str, doLoo: bool) -> Tuple[np.ndarray, np.ndarray]:
    x, A = computeAmplitude(dataFile, doLoo, field="rms")
    # print(x,A)
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
    w: float,
    x: np.ndarray,
    Nx: np.ndarray,
    x_flat: np.ndarray,
    Nx_flat: np.ndarray,
    x_start: float = 50,
    x_end: float = 250,
) -> float:
    """
    Computes the difference in N factor between a given case and a flat plate case at in a prescribed window of x values.
    """

    # interpolate x_flat, Nx_flat based on the values of x
    Nx_flat_interp = np.interp(x, x_flat, Nx_flat)

    # compute the average of deltaNx in the window w+50 < x < w+250
    # x_start = w + 80
    x_start = w + 80
    x_end = x_start + 50
    # x_start = findXStartMonotonicGrowth(w, x, Nx)
    # x_end = x_start + 200
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
        reynolds: np.ndarray | None = None,
        target_reynolds: float | None = None,
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
        self.reynolds = reynolds
        self.target_reynolds = target_reynolds


class SubsetKernel(Kernel):
    """Wraps a kernel so it only "sees" a chosen subset of input columns.

    sklearn's anisotropic/ARD kernels (e.g. ``Matern(length_scale=[...])``)
    give every input dimension its own length scale, but all dimensions are
    still combined into a *single* scaled Euclidean distance before applying
    one kernel function -- so every dimension is forced through the same
    kernel family/smoothness (nu). That's a reasonable assumption when all
    dimensions are "the same kind of thing" (width and depth both being
    spatial gap coordinates), but it's a much stronger assumption for
    Reynolds number: physically a global flow parameter, not a spatial
    coordinate, and in this dataset sampled at only 3 distinct values
    versus a dense sweep in (w, d).

    ``SubsetKernel`` is what lets us instead build a genuine
    product-of-subspaces kernel,

        k((w,d,Re),(w',d',Re')) = k_spatial(w,d ; w',d') * k_Re(Re ; Re'),

    so (w, d) and Re each get their own kernel family, length scale, and
    bounds, fit independently. get_params/theta/bounds/hyperparameters
    below mirror exactly how sklearn's own ``Sum``/``Product`` operators
    forward these to their children, so ``*``, cloning, and
    ``n_restarts_optimizer`` all work the same as for any built-in kernel.
    """

    def __init__(self, kernel: Kernel, dims):
        self.kernel = kernel
        self.dims = tuple(dims)

    def get_params(self, deep=True):
        params = dict(kernel=self.kernel, dims=self.dims)
        if deep:
            params.update(
                ("kernel__" + k, v) for k, v in self.kernel.get_params().items()
            )
        return params

    @property
    def hyperparameters(self):
        return [
            Hyperparameter("kernel__" + h.name, h.value_type, h.bounds, h.n_elements)
            for h in self.kernel.hyperparameters
        ]

    @property
    def theta(self):
        return self.kernel.theta

    @theta.setter
    def theta(self, theta):
        self.kernel.theta = theta

    @property
    def bounds(self):
        return self.kernel.bounds

    def __eq__(self, b):
        return type(self) == type(b) and self.kernel == b.kernel and self.dims == b.dims

    def _slice(self, X):
        return np.asarray(X)[:, self.dims]

    def __call__(self, X, Y=None, eval_gradient=False):
        Xs = self._slice(X)
        Ys = self._slice(Y) if Y is not None else None
        return self.kernel(Xs, Ys, eval_gradient=eval_gradient)

    def diag(self, X):
        return self.kernel.diag(self._slice(X))

    def is_stationary(self):
        return self.kernel.is_stationary()

    @property
    def requires_vector_input(self):
        return self.kernel.requires_vector_input

    def __repr__(self):
        return f"Subset(dims={self.dims}, {self.kernel})"


class ReynoldsAwareScaler:
    """Scales (width, depth) with an ordinary StandardScaler, and Reynolds
    with its own mean/std computed from the *distinct* Re levels present
    (not the raw, duplicate-heavy sample).

    With one dense campaign (Re=1000) and two sparse ones (Re=800, 3000),
    a plain StandardScaler fit on the raw Re column lets whichever Re has
    the most rows dominate the mean/std -- so "1 length-scale unit" in Re
    would depend on how many (w,d) points happened to be run at each Re,
    not on the physical spread of Re itself. Scaling from the unique
    levels fixes that. Exposes the same ``.transform(X)`` interface
    (X columns: width, depth, reynolds) as the plain StandardScaler it
    replaces.
    """

    def __init__(self):
        self.spatial_scaler = StandardScaler()
        self.re_mean_ = None
        self.re_std_ = None

    def fit(self, widths, depths, reynolds):
        self.spatial_scaler.fit(np.column_stack((widths, depths)))
        unique_re = np.unique(reynolds)
        self.re_mean_ = float(np.mean(unique_re))
        self.re_std_ = float(np.std(unique_re)) or 1.0
        return self

    def transform(self, X):
        X = np.asarray(X, dtype=float)
        sp = self.spatial_scaler.transform(X[:, :2])
        re = ((X[:, 2] - self.re_mean_) / self.re_std_).reshape(-1, 1)
        return np.column_stack((sp, re))


def _fit_gpr(
    depths: np.ndarray,
    widths: np.ndarray,
    reynolds: np.ndarray,
    deltaNx: np.ndarray,
    target_reynolds: float,
    reference_reynolds: float = 1000.0,
    w_max: float = 90.0,
    re_length_scale_bounds: Tuple[float, float] = (0.2, 30.0),
    alpha_real: float = 0.03,
    alpha_boundary: float = 0.01,
    alpha_prior: float = 0.6,
    add_prior: bool = True,
    n_boundary: int = 15,
    n_restarts_optimizer: int = 20,
    random_state: int = 42,
    verbose: bool = True,
):
    """
    Core GPR fit shared by ``gpr_interpolateWithRe`` and
    ``loo_cv_sparse_reynolds``. Builds a product kernel
    ``k_spatial(w,d) * k_Re(Re)`` (see ``SubsetKernel``) plus a
    ``WhiteKernel``, and uses a *per-point* alpha (noise variance) instead
    of one learned noise level shared by everything, so the near-exact
    zero-boundary points, the real CFD-derived points, and the crude
    empirical prior points don't compete on equal footing: boundary points
    are treated as essentially exact, real data gets a small
    measurement-noise-like variance, and the soft prior gets a
    deliberately large variance so it only weakly regularizes sparse
    regions instead of anchoring them (previously the prior points were
    commented as having "added noise" but none was actually applied).

    Returns ``(gp, scaler, orig_widths, orig_depths, orig_reynolds, orig_deltaNx)``.
    """
    depths = np.asarray(depths, dtype=float).ravel().copy()
    widths = np.asarray(widths, dtype=float).ravel().copy()
    reynolds = np.asarray(reynolds, dtype=float).ravel().copy()
    deltaNx = np.asarray(deltaNx, dtype=float).ravel().copy()

    if not (depths.size == widths.size == reynolds.size == deltaNx.size):
        raise ValueError(
            "depths, widths, reynolds, and deltaNx must have the same length"
        )

    valid_mask = (
        np.isfinite(depths)
        & np.isfinite(widths)
        & np.isfinite(reynolds)
        & np.isfinite(deltaNx)
    )
    depths, widths, reynolds, deltaNx = (
        depths[valid_mask],
        widths[valid_mask],
        reynolds[valid_mask],
        deltaNx[valid_mask],
    )
    if deltaNx.size == 0:
        raise ValueError("No valid DeltaN samples were provided")

    orig_widths, orig_depths = widths.copy(), depths.copy()
    orig_reynolds, orig_deltaNx = reynolds.copy(), deltaNx.copy()

    unique_reynolds = np.unique(reynolds)
    n_real = deltaNx.size

    # Zero-value boundary conditions (deltaN = 0 when w=0 or d=0) hold at
    # ANY Reynolds number -- "no gap" trivially means "no correction to N"
    # regardless of Re -- so it's valid physics to add them at every
    # observed Re plus the target Re, same as the original.
    max_depth = np.max(depths)
    boundary_reynolds = np.unique(np.append(unique_reynolds, target_reynolds))

    boundary_d = np.tile(np.linspace(0, max_depth, n_boundary), boundary_reynolds.size)
    boundary_w = np.zeros_like(boundary_d)
    boundary_re = np.repeat(boundary_reynolds, n_boundary)

    boundary_w2 = np.tile(np.linspace(0, w_max, n_boundary), boundary_reynolds.size)
    boundary_d2 = np.zeros_like(boundary_w2)
    boundary_re2 = np.repeat(boundary_reynolds, n_boundary)

    all_boundary_w = np.concatenate((boundary_w, boundary_w2))
    all_boundary_d = np.concatenate((boundary_d, boundary_d2))
    all_boundary_re = np.concatenate((boundary_re, boundary_re2))
    n_boundary_total = all_boundary_w.size

    widths_aug = np.concatenate((widths, all_boundary_w))
    depths_aug = np.concatenate((depths, all_boundary_d))
    reynolds_aug = np.concatenate((reynolds, all_boundary_re))
    deltaNx_aug = np.concatenate((deltaNx, np.zeros(n_boundary_total)))

    n_prior = 0
    if add_prior:
        # Soft prior points from an empirical model, anchored at
        # reference_reynolds only: this formula was calibrated against
        # (presumably) the dense Re=1000 data, so it isn't assumed to hold
        # at other Re -- that Re-dependence is exactly what's being fit.
        prior_w = np.array([60, 70, 80, 50, 60, 70, 80], dtype=float)
        prior_d = np.array([2.0, 2.0, 2.0, 2.5, 2.5, 2.5, 2.5])
        prior_deltaN = 0.1 * prior_w * np.tanh(44 * prior_d / prior_w)
        prior_reynolds = np.full_like(prior_w, reference_reynolds)
        widths_aug = np.concatenate((widths_aug, prior_w))
        depths_aug = np.concatenate((depths_aug, prior_d))
        reynolds_aug = np.concatenate((reynolds_aug, prior_reynolds))
        deltaNx_aug = np.concatenate((deltaNx_aug, prior_deltaN))
        n_prior = prior_w.size

    # Per-point noise variance: boundary << real data < soft prior.
    alpha = np.concatenate(
        [
            np.full(n_real, alpha_real**2),
            np.full(n_boundary_total, alpha_boundary**2),
            np.full(n_prior, alpha_prior**2),
        ]
    )

    scaler = ReynoldsAwareScaler().fit(widths_aug, depths_aug, reynolds_aug)
    X_train_scaled = scaler.transform(
        np.column_stack((widths_aug, depths_aug, reynolds_aug))
    )
    y_train = deltaNx_aug

    # Product kernel: (w,d) and Re get independent families/length-scales.
    k_spatial = SubsetKernel(
        Matern(length_scale=[1.0, 1.0], nu=2.5, length_scale_bounds=(1e-2, 1e2)),
        dims=[0, 1],
    )
    k_re = SubsetKernel(
        Matern(length_scale=1.0, nu=1.5, length_scale_bounds=re_length_scale_bounds),
        dims=[2],
    )
    kernel = k_spatial * k_re + WhiteKernel(
        noise_level=1e-2, noise_level_bounds=(1e-10, 1e-1)
    )

    gp = GaussianProcessRegressor(
        kernel=kernel,
        alpha=alpha,
        n_restarts_optimizer=n_restarts_optimizer,
        normalize_y=True,
        random_state=random_state,
    )
    gp.fit(X_train_scaled, y_train)

    if verbose:
        print(colors.OKBLUE + f"Learned kernel: {gp.kernel_}" + colors.ENDC)
        print(
            colors.OKBLUE
            + f"Log-marginal-likelihood: {gp.log_marginal_likelihood_value_:.3f}"
            + colors.ENDC
        )
        try:
            re_ls_scaled = gp.kernel_.k1.k2.kernel.length_scale
            print(
                colors.OKBLUE
                + f"Learned Re length-scale: {re_ls_scaled * scaler.re_std_:.0f} "
                + f"(Re units, from {unique_reynolds.size} distinct Re level(s))"
                + colors.ENDC
            )
        except AttributeError:
            pass

    return gp, scaler, orig_widths, orig_depths, orig_reynolds, orig_deltaNx


def gpr_interpolateWithRe(
    depths: np.ndarray,
    widths: np.ndarray,
    reynolds: np.ndarray,
    deltaNx: np.ndarray,
    target_reynolds: float | None = None,
    reference_reynolds: float = 1000.0,
    w_max: float = 90.0,
    grid_size: int = 200,
    re_length_scale_bounds: Tuple[float, float] = (0.2, 30.0),
    alpha_real: float = 0.03,
    alpha_boundary: float = 0.01,
    alpha_prior: float = 0.6,
    add_prior: bool = True,
    n_boundary: int = 15,
    n_restarts_optimizer: int = 20,
) -> GPRResult:
    """
    Fits a GPR model over (width, depth, Reynolds) and returns a 2D
    width-depth slice predicted at ``target_reynolds``.

    Unlike a single isotropic-per-dimension kernel over all three raw
    inputs, (w, d) and Re are given structurally different kernel
    treatment: a product kernel ``k_spatial(w,d) * k_Re(Re)`` (see
    ``SubsetKernel``), each with its own length scale and bounds, plus
    per-point noise that distinguishes real data, near-exact boundary
    conditions, and the soft empirical prior (see ``_fit_gpr``). This
    matters here specifically because Re is sampled at only a few distinct
    levels (dense at ``reference_reynolds``, sparse elsewhere), which a
    single combined-distance kernel doesn't treat any differently from a
    densely-swept spatial coordinate.

    If ``target_reynolds`` is not given, uses ``reference_reynolds`` when
    it is present in the training data, otherwise the median Reynolds
    value. ``re_length_scale_bounds`` are in units of "std dev of the
    *distinct* Re levels present" (see ``ReynoldsAwareScaler``) -- widen
    them if you have many more Re levels than the 3 this was tuned
    against, or if fits keep pinning the Re length-scale to a bound (check
    the printed "Learned Re length-scale" line).

    ``alpha_real``/``alpha_boundary``/``alpha_prior`` are noise *standard
    deviations* in DeltaN units. ``alpha_real`` in particular is a
    placeholder for your actual point-to-point scatter -- tighten or
    loosen it if you have a better estimate (e.g. from mesh-refinement or
    repeat runs).

    See ``loo_cv_sparse_reynolds`` to check how well this generalizes on
    your own data at the sparse Re levels, and
    ``plot_deltaN_vs_reynolds`` to inspect the learned Re-trend directly.
    """
    reynolds_arr = np.asarray(reynolds, dtype=float).ravel()
    unique_reynolds = np.unique(reynolds_arr[np.isfinite(reynolds_arr)])
    if target_reynolds is None:
        if np.any(np.isclose(unique_reynolds, reference_reynolds)):
            target_reynolds = reference_reynolds
        else:
            target_reynolds = float(np.median(unique_reynolds))
    else:
        target_reynolds = float(target_reynolds)

    print(
        colors.OKBLUE
        + f"Predicting DeltaN slice at Re_delta* = {target_reynolds:g}"
        + colors.ENDC
    )

    gp, scaler, orig_widths, orig_depths, orig_reynolds, orig_deltaNx = _fit_gpr(
        depths,
        widths,
        reynolds,
        deltaNx,
        target_reynolds=target_reynolds,
        reference_reynolds=reference_reynolds,
        w_max=w_max,
        re_length_scale_bounds=re_length_scale_bounds,
        alpha_real=alpha_real,
        alpha_boundary=alpha_boundary,
        alpha_prior=alpha_prior,
        add_prior=add_prior,
        n_boundary=n_boundary,
        n_restarts_optimizer=n_restarts_optimizer,
    )

    X = np.linspace(0, w_max, grid_size)
    Y = np.linspace(0, np.max(orig_depths), grid_size)
    X_grid, Y_grid = np.meshgrid(X, Y)
    XY = np.column_stack(
        (X_grid.ravel(), Y_grid.ravel(), np.full(X_grid.size, target_reynolds))
    )
    XY_scaled = scaler.transform(XY)

    Z_pred, sigma = gp.predict(XY_scaled, return_std=True)
    Z_grid = np.clip(Z_pred.reshape(X_grid.shape), 0, None)
    sigma_grid = sigma.reshape(X_grid.shape)

    return GPRResult(
        X_grid,
        Y_grid,
        Z_grid,
        sigma_grid,
        orig_widths,
        orig_depths,
        orig_deltaNx,
        gp,
        scaler,
        orig_reynolds,
        target_reynolds,
    )


def loo_cv_sparse_reynolds(
    depths: np.ndarray,
    widths: np.ndarray,
    reynolds: np.ndarray,
    deltaNx: np.ndarray,
    reference_reynolds: float = 1000.0,
    n_restarts_optimizer: int = 10,
    **fit_kwargs,
) -> list[dict]:
    """
    Leave-one-out cross-validation restricted to points *away* from
    ``reference_reynolds`` -- i.e. the sparse Re=800/3000-style campaigns.

    This is the part of the model the dense reference-Re data can't
    validate on its own, and the part most sensitive to how Re is treated
    in the kernel, so it's the most direct check of whether the
    Re-treatment actually generalizes on your real data (this was checked
    on synthetic data with a similar dense+sparse structure before this
    was handed back to you -- run this on your real data too).

    Each fold refits the whole model with one sparse point held out and
    predicts at that point, so this reuses the same n_restarts_optimizer
    per fold as a real fit; keep it modest (default 10) if you have many
    sparse points to loop over.

    Returns a list of dicts with keys: w, d, reynolds, actual, predicted,
    std, error.
    """
    depths = np.asarray(depths, dtype=float).ravel()
    widths = np.asarray(widths, dtype=float).ravel()
    reynolds = np.asarray(reynolds, dtype=float).ravel()
    deltaNx = np.asarray(deltaNx, dtype=float).ravel()

    sparse_idx = np.where(~np.isclose(reynolds, reference_reynolds))[0]
    if sparse_idx.size == 0:
        raise ValueError(
            "No points found away from reference_reynolds to cross-validate against"
        )

    fit_kwargs.setdefault("verbose", False)
    records = []
    keep = np.ones(len(widths), dtype=bool)
    for i in sparse_idx:
        keep[:] = True
        keep[i] = False

        gp, scaler, *_ = _fit_gpr(
            depths[keep],
            widths[keep],
            reynolds[keep],
            deltaNx[keep],
            target_reynolds=float(reynolds[i]),
            reference_reynolds=reference_reynolds,
            n_restarts_optimizer=n_restarts_optimizer,
            **fit_kwargs,
        )
        Xq = scaler.transform(np.array([[widths[i], depths[i], reynolds[i]]]))
        mean, std = gp.predict(Xq, return_std=True)
        records.append(
            dict(
                w=float(widths[i]),
                d=float(depths[i]),
                reynolds=float(reynolds[i]),
                actual=float(deltaNx[i]),
                predicted=float(mean[0]),
                std=float(std[0]),
                error=float(mean[0] - deltaNx[i]),
            )
        )

    errs = np.array([r["error"] for r in records])
    print(
        colors.OKBLUE
        + f"LOO-CV over {len(records)} sparse-Re point(s): "
        + f"MAE={np.mean(np.abs(errs)):.4f}, RMSE={np.sqrt(np.mean(errs**2)):.4f}, "
        + f"bias={np.mean(errs):+.4f}"
        + colors.ENDC
    )
    return records


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

    c = ax.pcolormesh(
        result.X_grid, result.Y_grid, result.Z_grid, shading="auto", cmap="viridis"
    )
    plot_mask = np.ones_like(result.deltaNx, dtype=bool)
    if result.reynolds is not None and result.target_reynolds is not None:
        plot_mask = np.isclose(result.reynolds, result.target_reynolds)

    ax.scatter(
        result.widths[plot_mask],
        result.depths[plot_mask],
        c=result.deltaNx[plot_mask],
        edgecolors="k",
        cmap="viridis",
        s=100,
    )

    for w, d, val in zip(
        result.widths[plot_mask], result.depths[plot_mask], result.deltaNx[plot_mask]
    ):
        ax.text(
            w,
            d - 0.05,
            f"{val:.2f}",
            color="white",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    contour_levels = [1, 2, 3, 4, 5]
    CS = ax.contour(
        result.X_grid,
        result.Y_grid,
        result.Z_grid,
        levels=contour_levels,
        colors="white",
        linewidths=1.2,
    )
    ax.clabel(CS, inline=True, fontsize=8, fmt="%d")

    fig.colorbar(c, ax=ax, label="Delta N")
    ax.set_xlabel("Width (w)")
    ax.set_ylabel("Depth (d)")
    title = "Interpolated and Extrapolated Delta N"
    if result.target_reynolds is not None:
        title += f" at Re_delta* = {result.target_reynolds:g}"
    ax.set_title(title)

    return fig


def plot_deltaN_vs_reynolds(
    result: GPRResult,
    points: list[tuple[float, float]],
    re_range: np.ndarray | None = None,
    ax: plt.Axes = None,
) -> plt.Figure:
    """
    Diagnostic plot: predicted DeltaN vs Reynolds at one or more fixed
    (w, d) locations, with a +/-2 sigma band -- the direct way to sanity
    check the learned Re-trend from ``gpr_interpolateWithRe`` (e.g. is it
    monotonic, does the band tighten near the Re you actually have data
    at, does it stay physically plausible when extrapolating).

    ``points`` is a list of (w, d) tuples. Any real data point within 5%
    of w_max in width and 5% of the max depth in depth of a requested
    (w, d) is overlaid as a marker at its own Reynolds number.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.get_figure()

    if re_range is None:
        re_lo = np.min(result.reynolds) if result.reynolds is not None else 800.0
        re_hi = np.max(result.reynolds) if result.reynolds is not None else 3000.0
        re_range = np.linspace(re_lo, re_hi, 60)

    w_tol = 0.05 * np.max(result.widths)
    d_tol = 0.05 * max(np.max(result.depths), 1.0)

    colors_cycle = plt.cm.tab10(np.linspace(0, 1, max(len(points), 1)))
    for (w0, d0), c in zip(points, colors_cycle):
        Xq = np.column_stack(
            (np.full_like(re_range, w0), np.full_like(re_range, d0), re_range)
        )
        Xq_scaled = result.scaler.transform(Xq)
        mean, std = result.gp.predict(Xq_scaled, return_std=True)
        mean = np.clip(mean, 0, None)

        ax.plot(re_range, mean, color=c, label=f"w={w0:g}, d={d0:g}")
        ax.fill_between(
            re_range,
            np.clip(mean - 2 * std, 0, None),
            mean + 2 * std,
            color=c,
            alpha=0.15,
        )

        if result.reynolds is not None:
            nearby = (np.abs(result.widths - w0) < w_tol) & (
                np.abs(result.depths - d0) < d_tol
            )
            if np.any(nearby):
                ax.scatter(
                    result.reynolds[nearby],
                    result.deltaNx[nearby],
                    color=c,
                    edgecolors="k",
                    zorder=5,
                    s=60,
                )

    ax.set_xlabel("Reynolds number")
    ax.set_ylabel("Delta N")
    ax.set_title("Predicted Delta N vs Re at fixed (w, d), +/-2 sigma")
    ax.legend()

    return fig
