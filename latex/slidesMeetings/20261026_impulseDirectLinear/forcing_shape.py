"""The forcing shape b(x, y) of the runs (session parameters x0 = -95, y0 = 0.5, weight w, Amplitude A = 0.5/w):

    b_rhou = -2 w (y - y0) A exp(-w r^2),   b_rhov = 2 w (x - x0) A exp(-w r^2),   r^2 = (x - x0)^2 + (y - y0)^2

i.e. b = (d psi/dy, -d psi/dx) with psi = A exp(-w r^2): a single Gaussian vortex (divergence free), counter-clockwise for A > 0.
"""
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID = "#1a1a19", "#5f5e5a", "#d9d8d2"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "savefig.dpi": 200})

x0, y0 = -95.0, 0.5
fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.9))
for ax, w in zip(axs, [1.0, 100.0]):
    A = 0.5 / w
    half = 3.0 / np.sqrt(w)
    x = np.linspace(x0 - half, x0 + half, 241)
    y = np.linspace(0.0, y0 + half, 161)  # fluid region y >= 0 (wall at y = 0)
    X, Y = np.meshgrid(x, y)
    G = A * np.exp(-w * ((X - x0) ** 2 + (Y - y0) ** 2))
    bu, bv = -2 * w * (Y - y0) * G, 2 * w * (X - x0) * G
    mag = np.hypot(bu, bv)
    cs = ax.contourf(X, Y, mag, levels=12, cmap="Blues")
    ax.streamplot(X, Y, bu, bv, color=INK, linewidth=0.7, density=0.9, arrowsize=0.7)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_title(f"$w$ = {w:g}: radius {1 / np.sqrt(2 * w):.2f}, max $|b|$ = {mag.max():.2f}",
                 fontsize=8, color=INK)
    fig.colorbar(cs, ax=ax, shrink=0.8, label="$|b|$")
axs[0].set_ylabel("$y$")
fig.tight_layout()
fig.savefig(f"{HERE}/Images/forcing_shape.png")
