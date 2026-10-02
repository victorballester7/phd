import numpy as np
from typing import Tuple
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter


def format_number(value: float) -> str:
    """Format run parameters without trailing .0 when they are integers."""
    if float(value).is_integer():
        return str(int(value))
    return f"{value:g}"


def timeFilter(
    time: np.ndarray, fields: np.ndarray, time_min: float, time_max: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Filters the time and fields arrays based on the specified time range.

    Args:
        time (np.ndarray): Array of time values.
    fields (np.ndarray): Array of field values (e.g., u, v, p) at every point (axis 0) and timestep (axis 1).
        time_min (float): Minimum time value to include in the comparison.
        time_max (float): Maximum time value to include in the comparison.
    Returns:
        tuple[np.ndarray, np.ndarray]: Filtered time and fields arrays.
    """

    # Ensure time is a 1D array
    if time.ndim != 1:
        raise ValueError("Time array must be 1-dimensional.")

    # Create a boolean mask for the time range
    mask = (time >= time_min) & (time <= time_max)

    # Filter the time and fields arrays using the mask
    filtered_time = time[mask]
    if fields.shape[0] == time.shape[0]:
        filtered_fields = fields[mask]
    else:
        filtered_fields = fields[:, mask]

    return filtered_time, filtered_fields


def spatialFilter(
    points: np.ndarray,
    fields: np.ndarray,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Filters the points and fields arrays based on the specified spatial range.

    Args:
        points (np.ndarray): Array of point coordinates (e.g., x, y, z) at every point (axis 0).
        fields (np.ndarray): Array of field values (e.g., u, v, p) at every point (axis 0) and timestep (axis 1).
        x_min (float): Minimum x-coordinate to include in the comparison.
        x_max (float): Maximum x-coordinate to include in the comparison.
        y_min (float): Minimum y-coordinate to include in the comparison.
        y_max (float): Maximum y-coordinate to include in the comparison.
    Returns:
        tuple[np.ndarray, np.ndarray]: Filtered points and fields arrays.
    """

    # Ensure points is a 2D array with at least 2 columns for x and y
    if points.ndim != 2 or points.shape[1] < 2:
        raise ValueError(
            "Points array must be 2-dimensional with at least 2 columns for x and y."
        )

    # Create a boolean mask for the spatial range
    mask = (
        (points[:, 0] >= x_min)
        & (points[:, 0] <= x_max)
        & (points[:, 1] >= y_min)
        & (points[:, 1] <= y_max)
    )

    # Filter the points and fields arrays using the mask
    filtered_points = points[mask]
    filtered_fields = fields[mask]

    return filtered_points, filtered_fields


def getQmean(data: np.ndarray) -> np.ndarray:
    """
    Mean kinetic-energy-like quantity Q = u^2 + v^2 averaged in time:

        Qbar = E[u^2 + v^2] = E[u]^2 + E[v]^2 + <u'u'> + <v'v'>

    where columns (u, v) of the avg file are the time-averaged velocities and
    (uu, uv, vv) are the Reynolds stresses (fluctuation covariances, mean
    already removed by Nektar's ReynoldsStresses filter).

    This is exactly ``getRMS(data)**2``, but it is the quantity that enters the
    L2 amplitude *linearly*, which is what makes the uncertainty propagation in
    ``pp.DeltaN_computation`` exact rather than linearized.
    """
    u = data[:, :, 2]
    v = data[:, :, 3]
    uu = data[:, :, 5]
    vv = data[:, :, 7]
    return np.abs(u**2 + v**2 + uu + vv)


def getRMS(data: np.ndarray) -> np.ndarray:
    """
    Computes the RMS of the fields array along the time axis, i.e. E(u^2 + v^2) = E(u)^2 + E(v)^2 + E(uu) + E(vv)
    """
    return np.sqrt(getQmean(data))


def getRMSVar(data: np.ndarray) -> np.ndarray:
    """
    Computes the variance of the RMS of the fields array along the time axis, i.e. Var(u^2 + v^2)

    Since the time series are gaussian and (u,v) are also jointly gaussian, using the Isserlis theorem we can compute the variance of the RMS as follows:
    Var(u^2 + v^2) = 4 E[u]^2E[uu] + 4 E[v]^2E[vv] + 2 E[uu]^2 + 2 E[vv]^2 + 8 E[u]E[v]E[uv] + 4 E[uv]^2

    Derivation (checked): write u = U + u', v = V + v' with (u', v') jointly
    gaussian, zero mean, covariances (<u'u'>, <u'v'>, <v'v'>) = (uu, uv, vv).
    Then Q = u^2 + v^2 = (U^2 + V^2) + (2U u' + 2V v') + (u'^2 + v'^2), so

        Var(Q) = Var(2U u' + 2V v') + Var(u'^2 + v'^2) + 2 Cov(linear, quadratic)

    The linear part gives 4 U^2 uu + 4 V^2 vv + 8 U V uv. The quadratic part
    gives, by Isserlis, Var(u'^2) = 2 uu^2, Var(v'^2) = 2 vv^2 and
    Cov(u'^2, v'^2) = 2 uv^2, i.e. 2 uu^2 + 2 vv^2 + 4 uv^2. The cross term
    vanishes because odd-order moments of a zero-mean gaussian are zero. Summing
    reproduces the formula above, so it is correct as written.

    IMPORTANT -- what this is and is not. This is the variance of the
    *instantaneous* signal Q(t) = u(t)^2 + v(t)^2, i.e. how much Q oscillates in
    time about its mean. It does NOT shrink as the averaging window grows, so on
    its own it cannot demonstrate convergence. The uncertainty of the *time
    average* that the .dat files actually store is

        Var(Qbar_estimated) = Var(Q) / n_eff,     n_eff = T_avg / (2 tau_int)

    with T_avg the averaging window and tau_int the integral time scale of Q.
    That one does decay like 1/T_avg. See ``effective_samples`` and the
    ``n_eff`` argument of ``pp.DeltaN_computation.computeAmplitude``.
    """
    u = data[:, :, 2]
    v = data[:, :, 3]
    uu = data[:, :, 5]
    uv = data[:, :, 6]
    vv = data[:, :, 7]
    var = (
        4 * u**2 * uu
        + 4 * v**2 * vv
        + 2 * uu**2
        + 2 * vv**2
        + 8 * u * v * uv
        + 4 * uv**2
    )
    return np.maximum(var, 0.0)


def effective_samples(T_avg: float, tau_int: float) -> float:
    """
    Effective number of independent samples in a time average of length
    ``T_avg`` for a stationary process with integral time scale ``tau_int``:

        n_eff = T_avg / (2 * tau_int)

    ``tau_int = int_0^inf rho(tau) dtau`` is the integral time scale of the
    *energy* signal Q = u^2 + v^2 (not of u itself), which for the narrowband
    TS-wave response to white-noise forcing is of the order of the inverse
    bandwidth of the amplified frequency band -- typically several TS periods,
    not one time step. Estimate it from the HistoryPoints time series of a run
    (autocorrelation of u^2 + v^2 at a point inside the boundary layer), then
    feed the resulting ``n_eff`` to the N-factor routines.

    Note that the raw sample count is a large overestimate of n_eff: with
    SampleFrequency=10 and dt=0.007 the samples are 0.07 apart, i.e. thousands
    of samples per TS period, and consecutive samples are almost perfectly
    correlated.
    """
    if tau_int <= 0:
        raise ValueError("tau_int must be positive")
    return float(T_avg) / (2.0 * float(tau_int))


def arc_length_parameterization(curve):
    """
    Arc-length parameter s in [0,1].
    """

    ds = np.sqrt(np.sum(np.diff(curve, axis=0) ** 2, axis=1))
    s = np.concatenate([[0.0], np.cumsum(ds)])
    s /= s[-1]

    return s


def average_curves(stable, unstable, npts=500):
    """
    Compute smooth stable/unstable boundaries by smoothing
    the middle curve and the gap separately.

    This guarantees that the two curves remain separated.
    """

    s_stable = arc_length_parameterization(stable)
    s_unstable = arc_length_parameterization(unstable)

    # ws_interp = PchipInterpolator(s_stable, stable[:, 0])
    # ds_interp = PchipInterpolator(s_stable, stable[:, 1])
    #
    # wu_interp = PchipInterpolator(s_unstable, unstable[:, 0])
    # du_interp = PchipInterpolator(s_unstable, unstable[:, 1])

    ws_interp = interp1d(s_stable, stable[:, 0], kind="linear")
    ds_interp = interp1d(s_stable, stable[:, 1], kind="linear")

    wu_interp = interp1d(s_unstable, unstable[:, 0], kind="linear")
    du_interp = interp1d(s_unstable, unstable[:, 1], kind="linear")

    # linear interpolation to avoid overshooting

    s = np.linspace(0.0, 1.0, npts)

    ws = ws_interp(s)
    ds = ds_interp(s)

    wu = wu_interp(s)
    du = du_interp(s)

    # ------------------------------------------------------
    # Middle curve
    # ------------------------------------------------------

    middle_x = 0.5 * (ws + wu)
    middle_y = 0.5 * (ds + du)
    middle_x = savgol_filter(middle_x, 41, 3)
    middle_y = savgol_filter(middle_y, 41, 3)

    # ------------------------------------------------------
    # Gap between branches
    # ------------------------------------------------------

    dgap = du - ds
    wgap = wu - ws

    # smooth the gap
    dgap_smooth = savgol_filter(
        dgap,
        window_length=21,  # must be odd
        polyorder=7,
    )

    wgap_smooth = savgol_filter(
        wgap,
        window_length=21,  # must be odd
        polyorder=7,
    )

    # guarantee positive separation
    # dgap_smooth = np.maximum(dgap_smooth, 1e-6)
    # wgap_smooth = np.maximum(wgap_smooth, 1e-6)

    # ------------------------------------------------------
    # Reconstruct branches
    # ------------------------------------------------------

    # plot wgap_smooth and dgap_smooth to check if they are positive

    # _,ax1 = plt.subplots()
    # _,ax2 = plt.subplots()
    # ax1.set_title(f"w_ini_stable={stable[0,0]}, d_ini_stable={stable[0,1]}")
    # ax2.set_title(f"w_ini_stable={stable[0,0]}, d_ini_stable={stable[0,1]}")
    # ax1.plot(s, du, label="du")
    # ax1.plot(s, ds, label="ds")
    # ax2.plot(s, wu, label="wu")
    # ax2.plot(s, ws, label="ws")

    stable_smooth = np.column_stack(
        (
            middle_x - 0.5 * wgap_smooth,
            middle_y - 0.5 * dgap_smooth,
        )
    )

    unstable_smooth = np.column_stack(
        (
            middle_x + 0.5 * wgap_smooth,
            middle_y + 0.5 * dgap_smooth,
        )
    )

    middle_curve = np.column_stack(
        (
            middle_x,
            middle_y,
        )
    )

    return stable_smooth, unstable_smooth, middle_curve
