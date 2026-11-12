"""Figures of the 20261026 meeting: compressible direct linear solver (Jacobian-free wrapper) and impulse response.

Data (interpolated with FieldConvert, see data/): pulse runs on the incompressible mesh (Ma 0.05, resolved dipole,
weight = 1) and impulse-response statistics (time integrals of u'^2) of the incompressible and compressible runs.
"""
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA, OUT = f"{HERE}/data", f"{HERE}/Images"

INK, INK2, GRID = "#1a1a19", "#5f5e5a", "#d9d8d2"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # validated categorical order (blue, orange, aqua, yellow)

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "legend.frameon": False, "lines.linewidth": 2, "savefig.dpi": 200,
})


def load(f):
    h = [c.strip() for c in open(f).readline().lstrip("#").split(",")]
    a = np.loadtxt(f, delimiter=",", skiprows=1)
    return {k: a[:, i] for i, k in enumerate(h)}


# ---------------------------------------------------------------------------------------------------------------
# 1. Pulse (wavepacket) on the incompressible mesh: before / after the rescaling-period fix
def pulse_figure():
    fig, axs = plt.subplots(2, 1, figsize=(7.2, 4.2), sharex=True)
    for ax, k in zip(axs, [3, 6]):
        inc = load(f"{DATA}/pulse_inc/chk_{k}.csv")
        old = load(f"{DATA}/pulse_commesh_old/chk_{k}.csv")
        new = load(f"{DATA}/pulse_verify/chk_{2 * k}.csv")
        m = inc["y"] > 0
        ax.plot(inc["x"][m], inc["u"][m], color=INK, lw=2.4, label="incompressible")
        ax.plot(new["x"][m], new["u"][m], color=SERIES[0], lw=1.6, ls="--",
                label=r"compressible, $T=24.5$, $\delta=10^{-6}$")
        ax.plot(old["x"][m], old["u"][m], color=SERIES[1], lw=1.4, ls=":",
                label=r"compressible, $T=98$, $\delta=10^{-5}$")
        ax.set_xlim(-100, 200)
        ax.set_ylabel(f"$u'$ at $y=0.5$\n$t={49 * k}$")
    axs[0].legend(loc="upper right", fontsize=8)
    axs[1].set_xlabel("$x$")
    axs[1].axvspan(0, 25, color=GRID, alpha=0.6, lw=0)
    axs[1].text(12.5, axs[1].get_ylim()[1] * 0.85, "gap", ha="center", color=INK2, fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{OUT}/pulse_fix.png")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------------
# 2. Impulse-response statistics: A(x) = max_y sqrt(int u'^2 dt), normalised by its value at x = -50
STATIONS = [-50, -20, 5, 12.5, 20, 30, 50, 100, 150, 200, 300, 400, 600, 800]


def integrals(f, T, inc):
    r = load(f)
    uu = (r["uu"] + (r["u"] ** 2 if inc else 0.0)) * T  # raw second moment x run time = int u'^2 dt
    return r["x"], r["y"], uu


def amplitude(f, T, inc):
    x, y, uu = integrals(f, T, inc)
    return np.array([np.sqrt(uu[(x == s) & (y != 0.5)].max()) for s in STATIONS])


def impulse_figure():
    I = f"{DATA}/impulse"
    cases = [
        ("incompressible", f"{I}/inc/mesh_37_stress.csv", 3626, True, dict(color=INK, lw=2.4, marker="o", ms=4)),
        ("Ma 0.05, incompressible mesh", f"{I}/incmesh/mesh_130_stats.csv", 3640, False,
         dict(color=SERIES[0], lw=1.8, ls="--", marker="s", ms=4)),
        ("Ma 0.05, DNS mesh", f"{I}/Ma0.05/mesh_130_stats.csv", 3640, False,
         dict(color=SERIES[1], lw=1.4, marker="^", ms=5)),
        ("Ma 0.2, DNS mesh", f"{I}/Ma0.2/mesh_130_stats.csv", 3640, False,
         dict(color=SERIES[2], lw=1.4, marker="v", ms=5)),
        ("Ma 0.4, DNS mesh", f"{I}/Ma0.4/mesh_130_stats.csv", 3640, False,
         dict(color=SERIES[3], lw=1.4, marker="D", ms=4)),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    table = {}
    for lab, f, T, inc, st in cases:
        A = amplitude(f, T, inc)
        table[lab] = A
        ax.semilogy(STATIONS, A / A[0], label=lab, **st)
    ax.axvspan(0, 25, color=GRID, alpha=0.6, lw=0)
    ax.set_xlabel("$x$")
    ax.set_ylabel(r"$\max_y \sqrt{\int u'^2 dt}$ / value at $x=-50$")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3)
    fig.tight_layout()
    fig.savefig(f"{OUT}/impulse_growth.png")
    plt.close(fig)
    print("normalised A(x):")
    print("  x     " + "".join(f"{k[:22]:>24s}" for k in table))
    for i, s in enumerate(STATIONS):
        print(f"{s:6.1f} " + "".join(f"{table[k][i] / table[k][0]:24.2f}" for k in table))


# 3. Wall-normal profiles of sqrt(int u'^2 dt), incompressible vs compressible on the same mesh
def profiles_figure():
    I = f"{DATA}/impulse"
    xi, yi, uui = integrals(f"{I}/inc/mesh_37_stress.csv", 3626, True)
    xc, yc, uuc = integrals(f"{I}/incmesh/mesh_130_stats.csv", 3640, False)
    fig, axs = plt.subplots(1, 4, figsize=(7.2, 2.8), sharey=True)
    for ax, s in zip(axs, [50, 150, 300, 600]):
        for x, y, uu, st, lab in [(xi, yi, uui, dict(color=INK, lw=2.4), "incompressible"),
                                  (xc, yc, uuc, dict(color=SERIES[0], lw=1.6, ls="--"), "Ma 0.05")]:
            m = (x == s) & (y != 0.5)
            p = np.sqrt(np.clip(uu[m], 0, None))
            ax.plot(p / p.max(), y[m], label=lab, **st)
        ax.set_title(f"$x={s}$", fontsize=9, color=INK)
        ax.set_xlabel(r"$\sqrt{\int u'^2}$ (norm.)", fontsize=8)
    axs[0].set_ylim(0, 8)
    axs[0].set_ylabel("$y$")
    axs[0].legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(f"{OUT}/impulse_profiles.png")
    plt.close(fig)


if __name__ == "__main__":
    pulse_figure()
    impulse_figure()
    profiles_figure()
