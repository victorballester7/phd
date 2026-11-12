"""What load_wave() (scripts/secondaryStabilityScan.py) reads, and what LocalWave does with it.

Diagnostic for docs/notes/20261108_secondaryStability.md.
Run from scripts/pp: .venv/bin/python plots/secondaryWaveInputs.py
Writes docs/notes/img/secondaryStability/waveInputs_w65-70.png and waveInputs_w60-64.png.
Rows: widths. Columns:
  (a) u' = u - U_BF   snapshot of the 2D limit cycle (2D run, pointsPO_n600.dat)
  (b) v'              same snapshot
  (c) U_BF            steady 2D base flow (baseflow2DwithW.fld on the same points)
  (d) at x - w = 150: the harmonics that LocalWave.at() fits to u', v' over two
      wavelengths, and the profiles that enter the stability operator
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../scripts"))
from secondaryStabilityScan import OMEGA, load_wave  # noqa: E402

OUT = os.path.join(HERE, "../../../docs/notes/img/secondaryStability")
X0 = 150


def figure(widths, out):
    fig, axs = plt.subplots(len(widths), 4, figsize=(19, 3.0 * len(widths)), squeeze=False,
                            gridspec_kw=dict(width_ratios=[1.3, 1.3, 1.3, 1]))
    for row, w in zip(axs, widths):
        L = load_wave(w)
        x, y = L.x, L.y
        jy = y <= 12
        ix = (x >= 0) & (x <= 450)
        for ax, q, name in ((row[0], L.u, "u' = u - U_BF"), (row[1], L.v, "v'"), (row[2], L.ub, "U_BF")):
            q = q[np.ix_(ix, jy)]
            if name == "U_BF":
                im = ax.pcolormesh(x[ix], y[jy], q.T, cmap="viridis", vmin=0, vmax=1, shading="auto")
                ax.contour(x[ix], y[jy], q.T, [0.99], colors="w", linewidths=0.6)
            else:
                m = np.abs(q).max()
                im = ax.pcolormesh(x[ix], y[jy], q.T, cmap="RdBu_r", vmin=-m, vmax=m, shading="auto")
            ax.axvline(X0, color="k", lw=0.6, ls="--")
            fig.colorbar(im, ax=ax, pad=0.01)
            ax.set_title(f"w={w}: {name}", fontsize=9)
            ax.set_ylabel("y")
        W = L.at(X0)
        ax = row[3]
        ax.plot(W["ubf"], y, "k", label="U_BF")
        ax.plot(W["ubf"] + W["u"][0].real, y, "k--", label="U_BF + u_0 (mean-flow distortion)")
        ax.axvline(W["c"], color="0.5", lw=0.6, ls=":", label=f"c = omega/alpha = {W['c']:.3f}")
        for n, col in ((1, "C3"), (2, "C1"), (3, "C2")):
            ax.plot(2 * np.abs(W["u"][n]), y, color=col, label=f"2|u_{n}|")
        ax.plot(20 * np.abs(W["v"][1]), y, "C0", label="20 x 2|v_1|")
        ax.set_ylim(0, 12)
        ax.set_xlim(0, 1.05)
        ax.set_title(f"w={w}, x-w={X0}: alpha={W['alpha']:.3f}, omega={OMEGA[w]}\n"
                     f"wave amplitude a = max 2|u_1| = {W['amp']:.3f}", fontsize=8)
        ax.legend(fontsize=6.5, loc="center right")
    for ax in axs[-1, :3]:
        ax.set_xlabel("x - w")
    axs[-1, 3].set_xlabel("velocity")
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print("saved", out)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    figure((65, 66, 70), f"{OUT}/waveInputs_w65-70.png")
    figure((60, 63, 64), f"{OUT}/waveInputs_w60-64.png")
