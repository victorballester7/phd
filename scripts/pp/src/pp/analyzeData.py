import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d


def shared_depth_curve(data_a, data_b):
    """
    Compute midpoint curve between two datasets at shared depth values.

    Uses vertical slices at common depth values to find the boundary
    between two stability regions.

    Parameters
    ----------
    data_a : np.ndarray
        First dataset (N×2 array with columns [w, d])
    data_b : np.ndarray
        Second dataset (N×2 array with columns [w, d])

    Returns
    -------
    tuple of np.ndarray
        (midpoint_w, midpoint_d) sorted by depth
    """
    # Extract unique depth values
    depths_a = np.unique(data_a[:, 1])
    depths_b = np.unique(data_b[:, 1])

    # Find intersection of depth values
    common_depths = np.intersect1d(depths_a, depths_b)

    midpoint_w = []
    midpoint_d = []

    for depth in common_depths:
        # Extract w values at this depth
        w_a = data_a[data_a[:, 1] == depth, 0]
        w_b = data_b[data_b[:, 1] == depth, 0]

        if len(w_a) == 0 or len(w_b) == 0:
            continue

        # Boundary points: rightmost of A, leftmost of B
        w_a_max = np.max(w_a)
        w_b_min = np.min(w_b)

        # Compute midpoint
        midpoint_w.append(0.5 * (w_a_max + w_b_min))
        midpoint_d.append(depth)

    # select max depth not in the common depths
    depth_above_all_stable = np.min(common_depths)
    width_max_depth_above_all_stable = np.min
    depth_all_stable = np.max(np.setdiff1d(depths_a, common_depths))

    # Convert to arrays and sort by depth
    midpoint_w = np.array(midpoint_w)
    midpoint_d = np.array(midpoint_d)
    sort_idx = np.argsort(midpoint_d)

    return midpoint_w[sort_idx], midpoint_d[sort_idx]


def build_bifurcation_pairs(data_set):
    """
    Build bifurcation pairs from a datasets.
    For example if data_set i made of 4 sub datasets as [eq, po1, po2, chaos],
    then we generate:
        (eq, po1), (po1, po2), (po2, chaos)
    """
    pairs = []

    for i in range(len(data_set) - 1):
        pairs.append((data_set[i], data_set[i + 1]))

    return pairs


# def compute_growthRate(time, data, ax, min_prominence=None, filter_outliers=True):
#     """
#     Compute growth rate of exponentially growing oscillations.

#     Uses envelope-based detrending and prominence filtering for robust
#     peak detection in signals with non-zero mean.

#     Parameters
#     ----------
#     time : np.ndarray
#         Array of time values (N,)
#     data : np.ndarray
#         Array of data values (N,)
#     ax : matplotlib axis
#         Axis for plotting (if needed)
#     min_prominence : float, optional
#         Minimum prominence for peak detection. If None, uses 2*std of detrended signal
#     filter_outliers : bool, optional
#         Whether to filter outliers using robust regression (default: True)

#     Returns
#     -------
#     tuple
#         (growth_rate, intercept) of the exponential envelope
#     """
#     
#     # Step 1: Estimate baseline/offset using median (robust to outliers)
#     baseline = np.median(data)
#     
#     # Step 2: Center data around baseline and take absolute value for envelope
#     data_centered = data - baseline
#     data_abs = np.abs(data_centered)
#     
#     # Step 3: Smooth absolute values to get approximate envelope (reduces noise)
#     window_size = max(5, len(data) // 50)  # adaptive window
#     if window_size % 2 == 0:
#         window_size += 1
#     from scipy.ndimage import uniform_filter1d
#     data_smoothed = uniform_filter1d(data_abs, size=window_size, mode='nearest')
#     
#     # Step 4: Find local maxima with prominence filtering
#     data_diff = np.diff(data_smoothed)
#     maxima_candidates = np.where((data_diff[:-1] > 0) & (data_diff[1:] < 0))[0] + 1
#     
#     # Remove boundary peaks
#     maxima_candidates = maxima_candidates[(maxima_candidates > 0) & 
#                                           (maxima_candidates < len(data) - 1)]
#     
#     if len(maxima_candidates) == 0:
#         return 0.0, np.log(np.mean(data_abs)) if np.mean(data_abs) > 0 else 0.0
#     
#     # Calculate prominence for each peak
#     prominences = []
#     valid_maxima = []
#     for idx in maxima_candidates:
#         # Find local minima on both sides
#         left_min = np.min(data_smoothed[max(0, idx-window_size):idx])
#         right_min = np.min(data_smoothed[idx+1:min(len(data_smoothed), idx+window_size+1)])
#         prominence = data_smoothed[idx] - max(left_min, right_min)
#         prominences.append(prominence)
#         valid_maxima.append(idx)
#     
#     prominences = np.array(prominences)
#     valid_maxima = np.array(valid_maxima)
#     
#     # Filter by prominence threshold
#     if min_prominence is None:
#         min_prominence = 2 * np.std(data_centered)
#     
#     significant_peaks = valid_maxima[prominences > min_prominence]
#     
#     if len(significant_peaks) < 2:
#         # Fallback: use top peaks if threshold too strict
#         n_peaks = min(len(valid_maxima), max(3, len(valid_maxima) // 3))
#         significant_peaks = valid_maxima[np.argsort(prominences)[-n_peaks:]]
#     
#     maxima = significant_peaks
#     
#     # Step 5: Interpolate to find precise peak times
#     times_max = time[maxima] - data_diff[maxima] * (
#         time[maxima + 1] - time[maxima]
#     ) / (data_diff[maxima + 1] - data_diff[maxima] + 1e-10)
#     
#     # Step 6: Get peak amplitudes from original (not smoothed) absolute data
#     data_max = np.interp(times_max, time, data_abs)
#     
#     # Filter out near-zero or negative values before log
#     valid_mask = data_max > 1e-10
#     if np.sum(valid_mask) < 2:
#         return 0.0, np.log(np.mean(data_abs[data_abs > 0])) if np.any(data_abs > 0) else 0.0
#     
#     times_max = times_max[valid_mask]
#     data_max = data_max[valid_mask]
#     
#     # Step 7: Log transform for exponential growth
#     log_data_max = np.log(data_max)
#     
#     # Step 8: Robust linear fitting with outlier rejection
#     if filter_outliers and len(times_max) > 3:
#         # First pass: standard fit
#         coeffs_initial = np.polyfit(times_max, log_data_max, 1)
#         predicted = np.polyval(coeffs_initial, times_max)
#         residuals = np.abs(log_data_max - predicted)
#         
#         # Filter outliers using median absolute deviation
#         mad = np.median(residuals)
#         threshold = 3 * mad
#         inliers = residuals < threshold
#         
#         if np.sum(inliers) >= 2:
#             times_max = times_max[inliers]
#             log_data_max = log_data_max[inliers]
#     
#     # Final fit with optional weighting (favor later peaks)
#     weights = np.linspace(0.5, 1.0, len(times_max))  # increasing weights
#     coefficients = np.polyfit(times_max, log_data_max, 1, w=weights)
#     growth_rate, intercept = coefficients
#     
#     return growth_rate, intercept

def compute_growthRate(time, data, ax):
    """
    Compute growth rate of a dataset over time.

    Parameters
    ----------
    time : np.ndarray
        Array of time values (N,)
    data : np.ndarray
        Array of data values (N,)

    Returns
    -------
    np.ndarray
        Array of growth rates (N-1,)
    """

    # data_abs = np.abs(data - np.mean(data))
    data_abs = np.abs(data - np.median(data))

    # get locations where data has a local maximum
    data_diff = np.diff(data_abs)
    maxima = np.where((data_diff[:-1] > 0) & (data_diff[1:] < 0))[0] + 1

    if maxima[-1] == len(data) - 1:
        maxima = maxima[:-1]

    # interpolation to find the time at which data reaches its maximum value
    times_max = time[maxima] - data_diff[maxima + 1] * (
        time[maxima + 1] - time[maxima]
    ) / (data_diff[maxima + 1] - data_diff[maxima])

    # interpolate data at times_max
    data_max = np.interp(times_max, time, data_abs)

    data_max = np.log(data_max)  # take log of data_max to get growth rate

    # line fitting to find the growth rate
    coefficients = np.polyfit(times_max, data_max, 1)
    growth_rate, intercept = coefficients

    return growth_rate, intercept
