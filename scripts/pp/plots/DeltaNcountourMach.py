"""
Delta n contour maps of the compressible gaps, one panel per Mach number (Re_delta* = 1000).

Data: the w = 1 impulse runs of the JacobianFreeDirect driver (R1 meshes, fixed base flows of the 2026-11-01
reruns), post-processed by scripts/compImpulseToPoints.py into
    src/comGapRe1000Ma<M>/directLinearSolver/blowingSuction/<dX_wY>/data/pointsimp_n600.dat
and, as the like-for-like reference, the flat-plate impulse run at the same Mach number (same forcing, same
numerics, uniform R1u mesh)
    src/flatPlateRe1000Ma<M>/directLinearSolver/blowingSuction/data/pointsimp_n600.dat

Delta n is computed exactly as for the incompressible map (pp.DeltaN_computation.computeDeltaN: mean of
N_gap - N_flat over x in [w+80, w+130]) and interpolated with the same GPR (gpr_interpolateWithRe, a single Re
level), but WITHOUT the empirical prior: its points (d = 2-2.5, w = 50-80) are an incompressible model and lie
in the unsteady region at every Mach number here. The zero boundary conditions at d = 0 and w = 0 are kept.

The hatched region is the 2D stability boundary of each Mach number (values from MaStabilityCurveManual.py),
extended at its last depth up to w = 90.

Outputs
    data/deltaN/DeltaNMa<M>_cases.dat   d, w, Delta n of every case
    data/deltaN/DeltaNMa<M>.dat         w, d, Delta n, sigma on the GPR grid
    images/DeltaNcountourMach.pdf       one panel per Mach number (+ the incompressible map if available)

usage: python plots/DeltaNcountourMach.py [--replot] [--mach 0.05 0.2 ...]
       --replot   reuse the .dat files instead of recomputing
"""

import argparse
import contextlib
import io
import os

import matplotlib.pyplot as plt
import numpy as np

from plots.deltaNstyle import COL_BOUNDARY, COL_EXCLUDED, deltaN_colormap_mach
from pp.colors import colors
from pp.DeltaN_computation import computeDeltaN, computeNx, gpr_interpolateWithRe
from pp.figureFrame import CM_TO_INCH, axis_on_top, jfm_ticks
from pp.fileManagement import extract_depth_width

plt.style.use("plots/style/jfm.mplstyle")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "../../..")
SRC = os.path.join(ROOT, "src")
DATA = os.path.join(ROOT, "data/deltaN")
FIG = os.path.join(ROOT, "images/DeltaNcountourMach.pdf")
N_PTS = 600
DATAFILE = f"pointsimp_n{N_PTS}.dat"
MACH = ["0.05", "0.2", "0.4", "0.6", "0.8"]
W_MAX, D_MAX = 90.0, 4.0
GRID = 200
CLIM = (0, 6)
LEVELS = [1, 2, 3, 4, 5]

# cases kept out of the fit (numerical variants of a case that is already in the list)
EXCLUDE = {"d0.5_w30_noSponge"}

# 2D stability boundary (d, w) per Mach number, copied from plots/MaStabilityCurveManual.py
BOUNDARY = {
    "0.05": [[4.0, 18.5], [3.5, 20.0], [3.0, 25.5], [2.5, 28.5], [2.0, 39.0], [1.5, 57.5]],
    "0.2": [[4.0, 16.5], [3.5, 17.5], [3.0, 18.5], [2.5, 24.5], [2.0, 29.0], [1.5, 37.5], [1.25, 50.0],
            [1.125, 60.0]],
    "0.4": [[4.0, 14.49], [3.5, 14.5], [3.0, 15.0], [2.5, 16.5], [2.0, 21.0], [1.5, 27.0], [1.0, 57.5],
            [0.875, 80.0]],
    "0.6": [[4.0, 12.5], [3.5, 13.49], [3.0, 13.5], [2.5, 14.5], [2.0, 19.0], [1.5, 27.5], [1.25, 35.5],
            [1.0, 52.5], [0.875, 60.0], [0.875, 80.0], [0.875, 120.0]],
    "0.8": [[4.0, 12.49], [3.5, 12.5], [3.0, 13.5], [2.5, 15.5], [2.0, 23.0], [1.5, 48.0], [1.25, 52.5],
            [1.125, 65.0], [1.125, 85.0]],
}


def boundary_wd(M):
    """Boundary as (w, d) sorted in w, extended at its last depth to W_MAX."""
    b = np.array(BOUNDARY[M], dtype=float)[:, ::-1]
    b = b[np.argsort(b[:, 0], kind="stable")]
    if b[-1, 0] < W_MAX:
        b = np.vstack((b, [W_MAX, b[-1, 1]]))
    return b[:, 0], b[:, 1]


def inside_unstable(M, w, d):
    bw, bd = boundary_wd(M)
    return w >= bw[0] and d > np.interp(w, bw, bd)


def quiet(f, *a):
    with contextlib.redirect_stdout(io.StringIO()):
        return f(*a)


def compute_mach(M):
    flat = os.path.join(SRC, f"flatPlateRe1000Ma{M}/directLinearSolver/blowingSuction/data", DATAFILE)
    if not os.path.isfile(flat):
        print(colors.WARNING + f"Ma {M}: no flat-plate reference ({flat}), skipped" + colors.ENDC)
        return None
    x_f, N_f = quiet(computeNx, flat, False)
    base = os.path.join(SRC, f"comGapRe1000Ma{M}/directLinearSolver/blowingSuction")
    rows = []
    for c in sorted(os.listdir(base)):
        f = os.path.join(base, c, "data", DATAFILE)
        if not os.path.isfile(f):
            continue
        d, w = extract_depth_width(c)
        x, N = quiet(computeNx, f, False)
        dN = quiet(computeDeltaN, w, x, N, x_f, N_f)
        tag = "excluded" if c in EXCLUDE else ("UNSTABLE SIDE" if inside_unstable(M, w, d) else "")
        print(f"  Ma {M:>4} {c:<20} Delta n = {dN:6.3f} {tag}")
        rows.append((d, w, dN, c in EXCLUDE))
    rows = np.array(rows, dtype=float)
    os.makedirs(DATA, exist_ok=True)
    np.savetxt(os.path.join(DATA, f"DeltaNMa{M}_cases.dat"), rows, fmt="%.6g",
               header="d w DeltaN excluded_from_fit  (impulse w=1, flat-plate impulse reference at the same Ma)")

    fit = rows[rows[:, 3] == 0]
    with contextlib.redirect_stdout(io.StringIO()):
        res = gpr_interpolateWithRe(fit[:, 0], fit[:, 1], np.full(len(fit), 1000.0), fit[:, 2],
                                    target_reynolds=1000.0, w_max=W_MAX, grid_size=GRID, add_prior=False)
    # predict on a common (w, d) grid for every Mach number (the GPR's own grid stops at the deepest case)
    W, D = np.meshgrid(np.linspace(0, W_MAX, GRID), np.linspace(0, D_MAX, GRID))
    XY = res.scaler.transform(np.column_stack((W.ravel(), D.ravel(), np.full(W.size, 1000.0))))
    Z, S = res.gp.predict(XY, return_std=True)
    Z = np.clip(Z.reshape(W.shape), 0, None)
    S = S.reshape(W.shape)
    np.savetxt(os.path.join(DATA, f"DeltaNMa{M}.dat"), np.column_stack((W.ravel(), D.ravel(), Z.ravel(), S.ravel())),
               fmt="%.6e", header=f"w d DeltaN sigma  ({GRID}x{GRID} grid, GPR without prior)")
    return W, D, Z, rows


def load_mach(M):
    g = os.path.join(DATA, f"DeltaNMa{M}.dat")
    c = os.path.join(DATA, f"DeltaNMa{M}_cases.dat")
    if not (os.path.isfile(g) and os.path.isfile(c)):
        return None
    a = np.loadtxt(g)
    n = int(round(np.sqrt(len(a))))
    return a[:, 0].reshape(n, n), a[:, 1].reshape(n, n), a[:, 2].reshape(n, n), np.atleast_2d(np.loadtxt(c))


def panel(ax, W, D, Z, title, M=None, rows=None, re_inc=False):
    c = ax.pcolormesh(W, D, Z, shading="auto", cmap=deltaN_colormap_mach(), edgecolors="none", linewidth=0,
                      antialiased=False, vmin=CLIM[0], vmax=CLIM[1], rasterized=True)
    cs = ax.contour(W, D, Z, levels=LEVELS, colors="black", linewidths=0.7, zorder=3)
    ax.clabel(cs, inline=True, fmt="%d", fontsize=7)
    ax.set_xlim(0, W_MAX)
    ax.set_ylim(0, D_MAX)
    if M is not None:
        bw, bd = boundary_wd(M)
        ax.fill_between(bw, bd, D_MAX, facecolor=COL_EXCLUDED, edgecolor="none", zorder=5)
        ax.plot(bw, bd, color=COL_BOUNDARY, ls=":", lw=0.9, zorder=6)
    if re_inc:
        from plots.deltaNstyle import add_excluded_region
        add_excluded_region(ax, 1000)
    if rows is not None:
        fit = rows[rows[:, 3] == 0]
        ax.plot(fit[:, 1], fit[:, 0], "o", ms=2.2, mfc="white", mec="black", mew=0.5, zorder=7)
    ax.set_title(title, pad=3)
    ax.grid(False)
    axis_on_top(ax)
    jfm_ticks(ax)
    return c


def plot(results):
    inc = os.path.join(DATA, "DeltaN1000.dat")
    panels = [(M, r) for M, r in results.items() if r is not None]
    n = len(panels) + (1 if os.path.isfile(inc) else 0)
    ncol = 3 if n > 4 else n
    nrow = int(np.ceil(n / ncol))
    w_cm, h_cm = 4.0, 3.2
    fig, axs = plt.subplots(nrow, ncol, figsize=((ncol * w_cm + 2.6) * CM_TO_INCH, (nrow * h_cm + 1.4) * CM_TO_INCH),
                            sharex=True, sharey=True, squeeze=False)
    fig.subplots_adjust(left=1.1 / (ncol * w_cm + 2.6), right=1 - 1.5 / (ncol * w_cm + 2.6),
                        bottom=0.9 / (nrow * h_cm + 1.4), top=1 - 0.5 / (nrow * h_cm + 1.4), wspace=0.12, hspace=0.25)
    flat = axs.ravel()
    letters = "abcdefgh"
    c = None
    for i, (M, (W, D, Z, rows)) in enumerate(panels):
        c = panel(flat[i], W, D, Z, rf"({letters[i]}) $\mathrm{{Ma}} = {M}$", M=M, rows=rows)
    if os.path.isfile(inc):
        a = np.loadtxt(inc)
        m = int(round(np.sqrt(len(a))))
        i = len(panels)
        c = panel(flat[i], a[:, 0].reshape(m, m), a[:, 1].reshape(m, m), a[:, 2].reshape(m, m),
                  rf"({letters[i]}) incompressible", re_inc=True)
    for ax in flat[n:]:
        ax.set_visible(False)
    for ax in axs[-1]:
        ax.set_xlabel(r"$w$")
    for ax in axs[:, 0]:
        ax.set_ylabel(r"$d$", rotation=0, labelpad=8)
        ax.set_yticks([0, 1, 2, 3, 4])
    cax = fig.add_axes([1 - 1.25 / (ncol * w_cm + 2.6), 0.9 / (nrow * h_cm + 1.4), 0.25 / (ncol * w_cm + 2.6),
                        1 - 1.4 / (nrow * h_cm + 1.4)])
    cb = fig.colorbar(c, cax=cax)
    cb.set_label(r"$\Delta n$", rotation=0, labelpad=8)
    fig.savefig(FIG, format="pdf")
    print(colors.OKGREEN + f"Figure saved to {FIG}" + colors.ENDC)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--replot", action="store_true")
    ap.add_argument("--mach", nargs="*", default=MACH)
    a = ap.parse_args()
    results = {M: (load_mach(M) if a.replot else compute_mach(M)) for M in a.mach}
    plot(results)


if __name__ == "__main__":
    main()
