"""
Sensitivity of the N factor curve of d4_w15 to the forcing used to excite the
perturbations (see ``scripts/makeNoiseForcing.py``):

    d4_w15           baseline: awgn, isotropic Gaussian, weight = 100
    d4_w15_envL      awgn, weight = 1
    d4_w15_envXL     awgn, weight = 0.2
    d4_w15_envStrip  awgn, weightX = 0.05, weightY = 4, x0 = -90
    d4_w15_msFlat    coherent multisine, flat PSD, w in [0.01, 0.25]
    d4_w15_msPink    coherent multisine, PSD ~ w^(-5/3)

All cases are averaged over the same window t in [5040, 24000]
(``combineAvg`` from the stress file at t >= 5000) and interpolated onto the
same ``data/pointsavg_n600.pts``. N(x) and DeltaN are computed exactly as in
``scripts/plotNfactorFromSections.py`` (L2 amplitude of u'^2 + v'^2, A0 the
minimum amplitude for x <= 0, DeltaN averaged over w+80 < x < w+130 against the
flat plate).

Top: N(x). Bottom: N(x) - N_baseline(x).
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from pp.DeltaN_computation import computeDeltaN, computeNx

plt.style.use("plots/style/jfm.mplstyle")

N = 600
DO_LOO = False
D, W = 4, 15

CASES = [
    # folder, label, colour, linestyle
    ("d4_w15", r"baseline: awgn, $w_g=100$", "black", "-"),
    ("d4_w15_envL", r"awgn, $w_g=1$", "tab:blue", "--"),
    ("d4_w15_envXL", r"awgn, $w_g=0.2$", "tab:cyan", "--"),
    ("d4_w15_envStrip", r"awgn, strip $0.05\times4$", "tab:green", "--"),
    ("d4_w15_msFlat", r"multisine, flat", "tab:orange", "-."),
    ("d4_w15_msPink", r"multisine, $\omega^{-5/3}$", "tab:red", "-."),
]


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    base = os.path.join(script_path, "../../../src/incGapRe1000/directLinearSolver/blowingSuction")
    flat_file = os.path.join(
        script_path, f"../../../src/flatPlateRe1000inc/directLinearSolver/blowingSuction/data/pointsavg_n{N}.dat"
    )

    x_flat, Nx_flat = computeNx(flat_file, DO_LOO)

    curves = {}
    for folder, label, c, ls in CASES:
        data_file = os.path.join(base, folder, "data", f"pointsavg_n{N}.dat")
        if not os.path.isfile(data_file):
            print(f"Missing {data_file}, skipping")
            continue
        x, Nx = computeNx(data_file, DO_LOO)
        dN = computeDeltaN(W, x, Nx, x_flat, Nx_flat)
        curves[folder] = (x, Nx, dN, label, c, ls)

    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(6.5, 6), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
    )

    ax.plot(x_flat, Nx_flat, color="gray", linewidth=3, alpha=0.5, label="flat plate")
    x_base, N_base = curves["d4_w15"][:2]
    for folder, (x, Nx, dN, label, c, ls) in curves.items():
        ax.plot(x, Nx, color=c, linestyle=ls, linewidth=1.2, label=rf"{label} ($\Delta N={dN:.2f}$)")
        if folder != "d4_w15":
            xc = x[(x >= x_base[0]) & (x <= x_base[-1])]
            ax2.plot(xc, np.interp(xc, x, Nx) - np.interp(xc, x_base, N_base), color=c, linestyle=ls, linewidth=1.2)

    for a in (ax, ax2):
        a.axvspan(0, W, color="0.85", zorder=0)  # gap
        a.axvspan(W + 80, W + 130, color="tab:purple", alpha=0.12, zorder=0)  # DeltaN window
    ax2.axhline(0, color="black", linewidth=0.8)

    ax.set_ylabel(r"$N$")
    ax2.set_ylabel(r"$N - N_{\mathrm{baseline}}$")
    ax2.set_xlabel(r"$x/\delta^*_0$")
    ax.set_title(rf"$d={D}$, $w={W}$: sensitivity to the forcing")
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()

    save_path = os.path.join(script_path, "../../../images/nfactorNoiseSource_d4_w15.pdf")
    fig.savefig(save_path, format="pdf")
    fig.savefig(save_path.replace(".pdf", ".png"), dpi=200)
    print(f"Saved {save_path}")

    print("\ncase               DeltaN    N(x=500)  N(x=1000)")
    for folder, (x, Nx, dN, *_) in curves.items():
        print(f"{folder:18s} {dN:7.3f}  {np.interp(500, x, Nx):8.3f}  {np.interp(1000, x, Nx):8.3f}")


if __name__ == "__main__":
    main()
