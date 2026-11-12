"""Leave-one-out check of the interpolation behind the per-Mach Delta n maps (plots/DeltaNcountourMach.py).

Reads data/deltaN/DeltaNMa<M>_cases.dat (d, w, Delta n, excluded) and compares, by leave-one-out prediction of every
case:
  A  per-Mach GPR        gpr_interpolateWithRe on one Mach number at a time (the current maps), no prior
  B  joint GPR           one GPR over (w, d, Ma): the existing third "Reynolds" input carries Ma, so every Mach
                         number borrows strength from its neighbours where it has no data
  C  model + GPR         per Mach, least-squares fit of the physical shape a w tanh(b d / w) (zero at d = 0 and
                         w = 0, linear in d for shallow gaps, saturating in w), then a per-Mach GPR on the residuals,
                         so away from the data the map falls back to the fitted shape instead of to a constant
Prints the LOO RMSE / max error per Mach number and overall, and the fitted (a, b) of model C.

usage: python scripts/deltaNMachGPRcheck.py [--data DIR] [--restarts 3]
"""

import argparse
import contextlib
import io
import os

import numpy as np
from scipy.optimize import curve_fit

from pp.DeltaN_computation import _fit_gpr

MACH = ["0.05", "0.2", "0.4", "0.6", "0.8"]
HERE = os.path.dirname(os.path.abspath(__file__))


def model(X, a, b):
    w, d = X
    return a * w * np.tanh(b * d / np.maximum(w, 1e-9))


def fit_model(w, d, y):
    p, _ = curve_fit(model, (w, d), y, p0=(0.05, 20.0), bounds=([0, 0.1], [1, 500]), maxfev=20000)
    return p


def gp_fit_predict(d, w, re, y, dq, wq, req, restarts):
    with contextlib.redirect_stdout(io.StringIO()):
        gp, sc, *_ = _fit_gpr(d, w, re, y, target_reynolds=float(np.median(req)), reference_reynolds=float(np.median(re)),
                              add_prior=False, n_restarts_optimizer=restarts, verbose=False)
    return gp.predict(sc.transform(np.column_stack((wq, dq, req))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "../../../data/deltaN"))
    ap.add_argument("--restarts", type=int, default=3)
    a = ap.parse_args()
    rows = []
    for M in MACH:
        f = os.path.join(a.data, f"DeltaNMa{M}_cases.dat")
        if os.path.isfile(f):
            r = np.atleast_2d(np.loadtxt(f))
            r = r[r[:, 3] == 0]
            rows += [(float(M), *x[:3]) for x in r]
    R = np.array(rows)  # Ma, d, w, dN
    ma, d, w, y = R.T
    err = {k: np.full(len(R), np.nan) for k in "ABC"}
    for i in range(len(R)):
        keep = np.arange(len(R)) != i
        s = keep & (ma == ma[i])
        # A
        err["A"][i] = gp_fit_predict(d[s], w[s], np.full(s.sum(), 1000.0), y[s], d[[i]], w[[i]], [1000.0],
                                     a.restarts)[0] - y[i]
        # B
        err["B"][i] = gp_fit_predict(d[keep], w[keep], ma[keep], y[keep], d[[i]], w[[i]], [ma[i]], a.restarts)[0] - y[i]
        # C
        p = fit_model(w[s], d[s], y[s])
        res = y[s] - model((w[s], d[s]), *p)
        err["C"][i] = (model((w[[i]], d[[i]]), *p) + gp_fit_predict(d[s], w[s], np.full(s.sum(), 1000.0), res, d[[i]],
                                                                       w[[i]], [1000.0], a.restarts))[0] - y[i]
        print(f"Ma {ma[i]:<4} d{d[i]:g}_w{w[i]:g}  dN {y[i]:6.3f}  LOO err A {err['A'][i]:+.3f} B {err['B'][i]:+.3f} "
              f"C {err['C'][i]:+.3f}", flush=True)
    print("\nLOO error          " + "".join(f"{k:>18}" for k in "ABC"))
    for M in [*np.unique(ma), None]:
        s = np.ones(len(R), bool) if M is None else ma == M
        lab = "all" if M is None else f"Ma {M:g} (n={s.sum()})"
        print(f"{lab:<18} " + "".join(f"   rmse {np.sqrt(np.mean(err[k][s] ** 2)):.3f} max {np.abs(err[k][s]).max():.2f}"
                                       for k in "ABC"))
    print("\nmodel C fits (all cases of each Mach): a w tanh(b d / w)")
    for M in np.unique(ma):
        s = ma == M
        p = fit_model(w[s], d[s], y[s])
        r = y[s] - model((w[s], d[s]), *p)
        print(f"  Ma {M:g}: a = {p[0]:.4f}, b = {p[1]:.2f}, residual rms {np.sqrt(np.mean(r ** 2)):.3f} "
              f"(incompressible model: a = 0.1, b = 44)")


if __name__ == "__main__":
    main()
