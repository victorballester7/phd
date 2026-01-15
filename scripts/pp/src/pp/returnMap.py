import numpy as np


def getReturnPoints(
        time: np.ndarray, field_section: np.ndarray, field_get_rp: np.ndarray, field_section_value: float, field_get_rp_min: float, field_get_rp_max: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Get the return points for a given section in phase space.

    Args:
        time (np.ndarray): Array of time values.
        field_section (np.ndarray): Array of field values to define the section (e.g., uL2).
        field_get_rp (np.ndarray): Array of field values from which to get the return points (e.g., vL2).

    Returns:
        tuple: A tuple containing:
        - return_times (np.ndarray): Array of time values at return points.
        - return_points (np.ndarray): Array of field values at return points.
    """
    print(f"Using section at field_section = {field_section_value} and getting return points from field_get_rp between {field_get_rp_min} and {field_get_rp_max}")
    
    # find indices where field_section crosses the section value
    tmp = field_section - field_section_value
    indices = np.where((tmp[:-1] * tmp[1:] < 0))[0] + 1

    # filter indices to keep only those where field_get_rp is within the specified range
    indices_field_rp_valid = np.where((field_get_rp[indices] >= field_get_rp_min) & (field_get_rp[indices] <= field_get_rp_max))[0]
    indices = indices[indices_field_rp_valid]
    
    # linear interpolation to find more accurate crossing points
    alpha = tmp[indices - 1] / (tmp[indices - 1] - tmp[indices])
    time_cross = time[indices - 1] + alpha * (time[indices] - time[indices - 1])
    field_cross = field_get_rp[indices - 1] + alpha * (field_get_rp[indices] - field_get_rp[indices - 1])

    return time_cross, field_cross

# # a: NxD array of modal amplitudes a(t) (D=2 or 3), times t (len N)
# # choose normal n (unit), level c, and direction sign_dir (1 or -1)

# s = a.dot(n) - c         # scalar distance from plane
# dsdt = np.gradient(s, t) # approximate derivative along trajectory

# crossings = []
# for k in range(len(s)-1):
#     if s[k]*s[k+1] < 0:                       # sign change => crossing
#         # enforce crossing direction:
#         if dsdt[k] * sign_dir <= 0:
#             continue
#         # linear interpolation to estimate tau where s = 0
#         alpha = s[k]/(s[k] - s[k+1])         # fraction from k to k+1
#         t_cross = t[k] + alpha*(t[k+1]-t[k])
#         a_cross = a[k] + alpha*(a[k+1]-a[k]) # interpolated point on section
#         # project a_cross onto a 1D coordinate along the section if needed:
#         # choose a basis vector u on the plane (orthonormal to n)
#         u = some_unit_vector_in_plane(n)
#         s_coord = u.dot(a_cross)
#         crossings.append((t_cross, a_cross, s_coord))

# # Build return map: s_n -> s_{n+1}
# s_coords = [c for (_,_,c) in crossings]
# # plot s_coords[:-1] vs s_coords[1:]
