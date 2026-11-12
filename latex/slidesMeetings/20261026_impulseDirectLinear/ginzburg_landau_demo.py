"""Impulse response vs white-noise forcing on the linearised Ginzburg-Landau equation (model of a convectively
unstable flow, Bagheri et al. 2009):

    dq/dt = L q + B xi(t),   L q = -nu q_x + gamma q_xx + mu(x) q,
    nu = U + 2 i c_u, gamma = 1 - i, mu(x) = mu0 - c_u^2 + mu2 x^2 / 2   (globally stable, convectively unstable)

Time stepping is exact for the linear part: q_{n+1} = A q_n + sqrt(dt) B w_n, A = exp(L dt), w_n i.i.d. N(0, 1)
(the same structure as the Nektar++ runs: a linear time-step map plus one noise sample per step).

Three ways to get the stationary covariance C = <q q^H>:
  1. exact:   discrete Lyapunov equation   C = A C A^H + dt B B^T
  2. impulse: C = dt sum_k (A^k B)(A^k B)^H  (one run from q_0 = B[:, j] per noise source j)
  3. noise:   time average of q q^H over a long forced run (statistical error ~ 1/sqrt(averaging time))

Cases: (a) one noise signal times a fixed shape b (rank 1): B = b;
       (b) independent noise at every grid point inside the Gaussian (like awgn() evaluated per quadrature point):
           B = diag(g), one source per point.
"""
import os

import matplotlib
import numpy as np
import scipy.linalg as sla

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID = "#1a1a19", "#5f5e5a", "#d9d8d2"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "legend.frameon": False, "lines.linewidth": 2, "savefig.dpi": 200,
})

# ----------------------------------------------------------------------------------------------- model
U, cu, mu0, mu2 = 2.0, 0.2, 0.3, -0.01  # max Re(lambda) = -0.075: globally stable
nu, gamma = U + 2j * cu, 1 - 1j
x = np.linspace(-60, 80, 281)
dx = x[1] - x[0]
n = len(x)
mu = mu0 - cu ** 2 + mu2 * x ** 2 / 2
D1 = (np.diag(np.ones(n - 1), 1) - np.diag(np.ones(n - 1), -1)) / (2 * dx)
D2 = (np.diag(np.ones(n - 1), 1) - 2 * np.eye(n) + np.diag(np.ones(n - 1), -1)) / dx ** 2
L = -nu * D1 + gamma * D2 + np.diag(mu)  # Dirichlet ends
lam = np.linalg.eigvals(L)
dt = 0.1
A = sla.expm(L * dt)

x0, sig = -11.0, 1.0
g = np.exp(-((x - x0) ** 2) / (2 * sig ** 2))  # forcing shape (Gaussian, like the dipole envelope)


def lyapunov(B):
    return sla.solve_discrete_lyapunov(A, dt * B @ B.conj().T)


def impulse(B, tol=1e-14):
    """dt * sum_k (A^k B)(A^k B)^H, accumulated only on the diagonal (pointwise <|q|^2>)."""
    var = np.zeros(n)
    Q = B.astype(complex).copy()
    q0 = np.sum(np.abs(Q) ** 2)
    k = 0
    while True:
        var += dt * np.sum(np.abs(Q) ** 2, axis=1)
        Q = A @ Q
        k += 1
        if np.sum(np.abs(Q) ** 2) < tol * q0 and k > 10:
            return var, k * dt
        if not np.isfinite(var).all() or k > 10 ** 6:
            raise RuntimeError("impulse response does not decay (globally unstable?)")


def noise_run(B, T_avg, T_spinup, seed):
    """Time average of |q|^2 over [T_spinup, T_spinup + T_avg] of the forced run; also the running averages."""
    rng = np.random.default_rng(seed)
    q = np.zeros(n, complex)
    m = B.shape[1]
    nsp, navg = int(T_spinup / dt), int(T_avg / dt)
    sq = np.sqrt(dt)
    for _ in range(nsp):
        q = A @ q + sq * (B @ rng.standard_normal(m))
    acc = np.zeros(n)
    checkpoints, running = [], []
    marks = set(np.unique(np.logspace(1, np.log10(navg), 25).astype(int)))
    for i in range(1, navg + 1):
        q = A @ q + sq * (B @ rng.standard_normal(m))
        acc += np.abs(q) ** 2
        if i in marks:
            checkpoints.append(i * dt)
            running.append(acc / i)
    return acc / navg, np.array(checkpoints), running


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


if __name__ == "__main__":
    assert lam.real.max() < 0, "globally unstable model"
    print(f"max Re(lambda) = {lam.real.max():.4f} (< 0: globally stable)")
    # (a) rank-1 forcing: one noise signal times the shape g
    Ba = g[:, None]
    Ca = np.clip(np.real(np.diag(lyapunov(Ba))), 0, None)
    Ia, Timp = impulse(Ba)
    print(f"(a) rank 1: impulse vs Lyapunov {rel(Ia, Ca):.2e}  (impulse run length {Timp:.0f})")
    T_spin = 300
    Na, tc, runs = noise_run(Ba, 20000, T_spin, seed=1)
    errs = [rel(r, Ca) for r in runs]
    print(f"(a) rank 1: noise run (T_avg = 20000) vs Lyapunov {rel(Na, Ca):.2e}")
    # (b) independent noise at each grid point (awgn per quadrature point)
    pts = np.where(g > 1e-3)[0]
    Bb = np.diag(g)[:, pts]
    Cb = np.clip(np.real(np.diag(lyapunov(Bb))), 0, None)
    Ib_all, _ = impulse(Bb)  # one impulse per source (column)
    Nb, _, _ = noise_run(Bb, 20000, T_spin, seed=2)
    print(f"(b) per point ({len(pts)} sources): sum of {len(pts)} impulses vs Lyapunov {rel(Ib_all, Cb):.2e}, "
          f"noise run vs Lyapunov {rel(Nb, Cb):.2e}, single impulse from g vs Lyapunov {rel(Ia, Cb):.2e}")
    s_a = np.sqrt(Ca / Ca.max())
    s_b = np.sqrt(Cb / Cb.max())
    dn = x > x0 + 5  # downstream of the forcing
    print(f"    shape difference (normalised rms, x > x0 + 5) rank-1 vs per-point: {rel(s_a[dn] / s_a[dn].max(), s_b[dn] / s_b[dn].max()):.2e}; level ratio "
          f"max C_a / max C_b = {Ca.max() / Cb.max():.2f}")

    # figure 1: rank-1 case, the three methods + convergence of the noise run
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.5))
    ax = axs[0]
    ax.semilogy(x, np.sqrt(Ca), color=INK, lw=2.6, label="exact (Lyapunov)")
    ax.semilogy(x, np.sqrt(Ia), color=SERIES[0], lw=1.6, ls="--", label=f"impulse (1 run, $t\\leq${Timp:.0f})")
    i1 = np.argmin(np.abs(tc - 1000))
    ax.semilogy(x, np.sqrt(runs[i1]), color=SERIES[1], lw=1.0, label="noise, $T_{avg}=10^3$")
    ax.semilogy(x, np.sqrt(Na), color=SERIES[2], lw=1.2, ls=":", label="noise, $T_{avg}=2\\times10^4$")
    ax.axvline(x0, color=INK2, lw=0.8, ls="--")
    ax.set_xlabel("$x$")
    ax.set_ylabel(r"$\sqrt{\langle|q'|^2\rangle}$")
    ax.set_ylim(np.sqrt(Ca.max()) * 1e-3, np.sqrt(Ca.max()) * 3)
    ax.legend(fontsize=7, loc="upper right")
    ax = axs[1]
    ax.loglog(tc, errs, color=SERIES[1], lw=1.6, marker="o", ms=3, label="noise run")
    ax.loglog(tc, errs[-1] * np.sqrt(tc[-1] / tc), color=INK2, lw=1, ls="--", label=r"$\propto T_{avg}^{-1/2}$")
    ax.axhline(rel(Ia, Ca), color=SERIES[0], lw=1.6, ls="--", label="impulse")
    ax.set_xlabel("averaging time $T_{avg}$")
    ax.set_ylabel("relative error vs exact")
    ax.legend(fontsize=7, loc="center right")
    fig.tight_layout()
    fig.savefig(f"{HERE}/Images/gl_rank1.png")
    plt.close(fig)

    # figure 2: per-point noise
    fig, ax = plt.subplots(figsize=(7.4, 2.3))
    ax.semilogy(x, np.sqrt(Cb), color=INK, lw=2.6, label="per-point noise: exact (Lyapunov)")
    ax.semilogy(x, np.sqrt(Nb), color=SERIES[2], lw=1.2, ls=":", label="per-point noise: noise run, $T_{avg}=2\\times10^4$")
    ax.semilogy(x, np.sqrt(Ib_all), color=SERIES[0], lw=1.6, ls="--", label=f"sum of {len(pts)} impulses (one per point)")
    ax.semilogy(x, np.sqrt(Ia), color=SERIES[1], lw=1.4, label="single impulse from the shape $g$ (= rank-1 forcing)")
    ax.axvline(x0, color=INK2, lw=0.8, ls="--")
    ax.set_xlabel("$x$")
    ax.set_ylabel(r"$\sqrt{\langle|q'|^2\rangle}$")
    ax.set_ylim(np.sqrt(Cb.max()) * 1e-3, np.sqrt(Ia.max()) * 3)
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(f"{HERE}/Images/gl_perpoint.png")
    plt.close(fig)
