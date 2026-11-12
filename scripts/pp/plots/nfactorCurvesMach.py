"""
N factor curves N(x) of the compressible gaps, one panel per Mach number (Re_delta* = 1000).

Same data as plots/DeltaNcountourMach.py: the w = 1 impulse runs post-processed by
scripts/compImpulseToPoints.py into
    src/comGapRe1000Ma<M>/directLinearSolver/blowingSuction/<dX_wY>/data/pointsimp_n600.dat
and the flat-plate impulse reference at the same Mach number
    src/flatPlateRe1000Ma<M>/directLinearSolver/blowingSuction/data/pointsimp_n600.dat

Every gap curve is coloured by its Delta n (same colour map and limits as the contour maps); the flat plate is
the dashed black curve. Numerical variants (EXCLUDE) are left out, and so is Ma0.05 d1.5_w25 (EXCLUDE_MACH): its
production impulse ran on the coarse DNS mesh (the 2026-10-16 pilot), not on R1, so the TS waves are not resolved
downstream and its N falls below the flat plate from x ~ 150 on (found 2026-11-12). The curves stop at X_MAX,
before the outlet region, where every N drops (boundary effect, from x = 935).

Output: images/nfactorCurvesMach.pdf
usage:  python plots/nfactorCurvesMach.py
"""

import contextlib
import io
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from pp.colors import colors
from pp.DeltaN_computation import computeDeltaN, computeNx
from pp.figureFrame import CM_TO_INCH, jfm_ticks
from pp.fileManagement import extract_depth_width

plt.style.use("plots/style/jfm.mplstyle")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "../../..")
SRC = os.path.join(ROOT, "src")
FIG = os.path.join(ROOT, "images/nfactorCurvesMach.pdf")
DATAFILE = "pointsimp_n600.dat"
MACH = ["0.05", "0.2", "0.4", "0.6", "0.8"]
EXCLUDE = {"d0.5_w30_noSponge"}
EXCLUDE_MACH = {("0.05", "d1.5_w25")}  # impulse on the coarse DNS mesh, not R1 (see above)
CLIM = (0, 5)
CMAP = mpl.colors.LinearSegmentedColormap.from_list(  # Oranges without its near-white end, so that low Delta n shows
    "OrangesLines", plt.get_cmap("Oranges")(np.linspace(0.25, 1.0, 256)))
X_MAX = 925


def nx(path):
    with contextlib.redirect_stdout(io.StringIO()):
        x, N = computeNx(path, False)
    keep = x <= X_MAX
    return x[keep], N[keep]


def load(M):
    flat = nx(os.path.join(SRC, f"flatPlateRe1000Ma{M}/directLinearSolver/blowingSuction/data", DATAFILE))
    base = os.path.join(SRC, f"comGapRe1000Ma{M}/directLinearSolver/blowingSuction")
    gaps = []
    for c in sorted(os.listdir(base)):
        f = os.path.join(base, c, "data", DATAFILE)
        if c in EXCLUDE or (M, c) in EXCLUDE_MACH or not os.path.isfile(f):
            continue
        d, w = extract_depth_width(c)
        x, N = nx(f)
        with contextlib.redirect_stdout(io.StringIO()):
            dN = computeDeltaN(w, x, N, *flat)
        gaps.append((dN, c, x, N))
    return flat, sorted(gaps)  # low Delta n first, so the strongest curves are drawn on top


def main():
    cmap = CMAP
    norm = mpl.colors.Normalize(*CLIM)
    ncol, nrow = 3, 2
    w_cm, h_cm = 3.15, 2.8  # total width 13.35 cm <= JFM textwidth 13.49 cm
    pad_l = pad_r = 1.6  # symmetric: the colorbar lives in pad_r
    pad_b, pad_t, gap_x, gap_y = 0.9, 0.5, 0.35, 0.65
    fig_w = pad_l + ncol * w_cm + (ncol - 1) * gap_x + pad_r
    fig_h = pad_b + nrow * h_cm + (nrow - 1) * gap_y + pad_t
    fig = plt.figure(figsize=(fig_w * CM_TO_INCH, fig_h * CM_TO_INCH))

    def axes_cm(left, bottom, width, height, **kw):
        return fig.add_axes([left / fig_w, bottom / fig_h, width / fig_w, height / fig_h], **kw)

    axs = []
    for k in range(ncol * nrow):
        i, j = divmod(k, ncol)
        axs.append(axes_cm(pad_l + j * (w_cm + gap_x), pad_b + (nrow - 1 - i) * (h_cm + gap_y), w_cm, h_cm,
                           sharex=axs[0] if axs else None, sharey=axs[0] if axs else None))

    letters = "abcdef"
    for k, M in enumerate(MACH):
        ax = axs[k]
        flat, gaps = load(M)
        for dN, c, x, N in gaps:
            ax.plot(x, N, color=cmap(norm(dN)), linewidth=0.8)
        ax.plot(*flat, color="black", linestyle="--", linewidth=1.2)
        ax.set_title(rf"({letters[k]}) $\mathrm{{Ma}} = {M}$ ({len(gaps)} gaps)", fontsize=8)
        print(f"Ma {M}: {len(gaps)} gaps, Delta n {min(g[0] for g in gaps):.2f}-{max(g[0] for g in gaps):.2f}, "
              f"N_flat({X_MAX}) = {flat[1][-1]:.2f}, max N_gap = {max(g[3].max() for g in gaps):.2f}")

    # last slot: key
    key = axs[-1]
    key.axis("off")
    key.plot([], [], color="black", linestyle="--", linewidth=1.2, label="flat plate")
    key.plot([], [], color=cmap(norm(3)), linewidth=0.8, label=r"gaps (colour: $\Delta n$)")
    key.legend(loc="center", frameon=False, fontsize=8)

    for k, ax in enumerate(axs[:-1]):
        ax.set_xlim(-70, X_MAX)
        jfm_ticks(ax)
        if k % ncol == 0:
            ax.set_ylabel(r"$n$", rotation=0, ha="right", va="center", labelpad=6)
        else:
            ax.tick_params(labelleft=False)
        if k >= len(MACH) - ncol:  # bottom panel of its column
            ax.set_xlabel(r"$x$")
        else:
            ax.tick_params(labelbottom=False)

    cax = axes_cm(fig_w - pad_r + 0.3, pad_b, 0.25, nrow * h_cm + (nrow - 1) * gap_y)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label(r"$\Delta n$", rotation=0, ha="left", va="center")
    fig.savefig(FIG, format="pdf")
    print(colors.OKGREEN + f"Figure saved to {FIG}" + colors.ENDC)


if __name__ == "__main__":
    main()
