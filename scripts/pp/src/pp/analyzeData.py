import numpy as np

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

