import numpy as np
from scipy.interpolate import make_smoothing_spline

def smooth_curve(w, d, lamb = 0.001):
    spline = make_smoothing_spline(d, w, lam=lamb)
    d_smooth = np.linspace(np.min(d), np.max(d), 500)
    w_smooth = spline(d_smooth)
    return w_smooth, d_smooth
