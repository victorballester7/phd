"""
Local secondary (Floquet) stability of a finite-amplitude 2D travelling wave,
in the spirit of Herbert (1988), with the wave taken from a DNS snapshot.

Base flow, in the frame moving with the wave at its phase speed c, parallel and
periodic in xi = x - c t with wavenumber alpha:

    U(xi, y) = sum_n U_n(y) e^{i n alpha xi},   V(xi, y) = sum_n V_n(y) e^{i n alpha xi},

where U_0 is the steady base flow plus the mean-flow distortion minus c, and V_0 = 0
(parallel-flow approximation). The 3D disturbance is

    q'(xi, y, z, t) = e^{sigma t + i beta z} sum_{m=-M}^{M} q_m(y) e^{i (m + gamma) alpha xi},

with gamma = 0 for fundamental (K-type) and gamma = 1/2 for subharmonic (H-type)
modes. The linearised Navier-Stokes equations in primitive variables (u, v, w, p)
are collocated on Chebyshev nodes in y in [0, ymax] (algebraic map clustering
half of the nodes below yi), with no-slip at the wall, u = v = w = 0 at ymax and
continuity at every node; the result is a generalised eigenvalue problem for the
temporal growth rate sigma (frame independent). The spatial rate follows as
sigma_x = Re(sigma) / c, since the fundamental mode travels with the wave.

Lengths in delta*_le, velocities in U_inf, Re = Re_delta*_le, as in the DNS.

Checks (scratch tests, 2026-11-05): with no wave, the harmonic m = 0, gamma = 1,
beta = 0 is the Orr-Sommerfeld problem and gives the Blasius critical point
(Re_delta* = 519.4, alpha = 0.3034: neutral, c = 0.3965); with a finite wave the
eigenvalues agree with an independent version collocated in xi; and they change
by < 0.5 % from (N, M, ymax) = (48, 3, 40) to (64, 5, 60).
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as sla
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
from scipy.signal import hilbert


def cheb(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Chebyshev-Gauss-Lobatto nodes s_j = cos(pi j / n) and differentiation matrix."""
    s = np.cos(np.pi * np.arange(n + 1) / n)
    c = np.r_[2.0, np.ones(n - 1), 2.0] * (-1) ** np.arange(n + 1)
    dS = s[:, None] - s[None, :]
    D = np.outer(c, 1 / c) / (dS + np.eye(n + 1))
    D -= np.diag(D.sum(1))
    return s, D


def grid(n: int, ymax: float, yi: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Mapped nodes y = a (1 + s) / (b - s) on [0, ymax] (node 0 at ymax), D and D^2."""
    s, Ds = cheb(n)
    a = yi * ymax / (ymax - 2 * yi)
    b = 1 + 2 * a / ymax
    y = a * (1 + s) / (b - s)
    D1 = np.diag((b - s) ** 2 / (a * (1 + b))) @ Ds
    return y, D1, D1 @ D1


class SecondaryStability:
    """Linearised operator about the wave; harmonics m = -M..M of the disturbance."""

    def __init__(self, y, D1, D2, Un: dict, Vn: dict, alpha: float, Re: float, M: int):
        self.y, self.D1, self.D2 = y, D1, D2
        self.Un, self.Vn, self.alpha, self.nu, self.M = Un, Vn, alpha, 1.0 / Re, M
        self.dUn = {n: D1 @ u for n, u in Un.items()}
        self.dVn = {n: D1 @ v for n, v in Vn.items()}

    def matrices(self, beta: float, gamma: float) -> tuple[np.ndarray, np.ndarray]:
        y, D1, D2, M, al, nu = self.y, self.D1, self.D2, self.M, self.alpha, self.nu
        ny = len(y)
        I = np.eye(ny)
        ms = np.arange(-M, M + 1)
        nb = 4 * ny
        A = np.zeros((len(ms) * nb, len(ms) * nb), complex)
        B = np.zeros(len(ms) * nb)
        k = (ms + gamma) * al

        def blk(im, var):  # var: 0 u, 1 v, 2 w, 3 p
            o = im * nb + var * ny
            return slice(o, o + ny)

        for im in range(len(ms)):
            lap = nu * (D2 - (k[im] ** 2 + beta ** 2) * I)
            ru, rv, rw, rc = blk(im, 0), blk(im, 1), blk(im, 2), blk(im, 3)
            # sigma q = -L q on the momentum rows
            A[ru, blk(im, 0)] += lap
            A[rv, blk(im, 1)] += lap
            A[rw, blk(im, 2)] += lap
            A[ru, blk(im, 3)] -= 1j * k[im] * I
            A[rv, blk(im, 3)] -= D1
            A[rw, blk(im, 3)] -= 1j * beta * I
            for n in self.Un:  # base harmonic n acting on disturbance harmonic m - n
                jm = im - n
                if not 0 <= jm < len(ms):
                    continue
                U, V, dU, dV = self.Un[n], self.Vn[n], self.dUn[n], self.dVn[n]
                adv = np.diag(1j * k[jm] * U) + np.diag(V) @ D1
                A[ru, blk(jm, 0)] -= adv + np.diag(1j * n * al * U)
                A[ru, blk(jm, 1)] -= np.diag(dU)
                A[rv, blk(jm, 1)] -= adv + np.diag(dV)
                A[rv, blk(jm, 0)] -= np.diag(1j * n * al * V)
                A[rw, blk(jm, 2)] -= adv
            A[rc, blk(im, 0)] = 1j * k[im] * I
            A[rc, blk(im, 1)] = D1
            A[rc, blk(im, 2)] = 1j * beta * I
            for var in range(3):
                r = blk(im, var)
                B[r] = 1.0
                for j in (r.start, r.stop - 1):  # Dirichlet at ymax and at the wall
                    A[j, :] = 0.0
                    A[j, j] = 1.0
                    B[j] = 0.0
        return A, B

    def eigenvalues(self, beta: float, gamma: float, nev: int = 1) -> np.ndarray:
        """Leading temporal eigenvalues sigma, sorted by decreasing growth rate."""
        A, B = self.matrices(beta, gamma)
        w = sla.eig(A, np.diag(B), right=False)
        w = w[np.isfinite(w) & (np.abs(w) < 1e3)]
        return w[np.argsort(-w.real)][:nev]


class LocalWave:
    """Harmonics of a 2D wave from a snapshot on a structured (x, y) set of points.

    u, v: perturbation (u - U_BF) on (x, y); ub: steady base flow U_BF on the same
    points; omega: frequency of the wave (from history points).
    """

    def __init__(self, x, y, u, v, ub, omega: float):
        self.x, self.y, self.u, self.v, self.ub, self.omega = x, y, u, v, ub, omega
        # local wavenumber from the phase of v along x at the height of its maximum
        j = np.argmax(np.abs(v).max(0) * (y < 10))
        ph = np.unwrap(np.angle(hilbert(v[:, j])))
        self.alpha_x = np.abs(np.gradient(ph, x))

    def _fit(self, x0: float, alpha: float, nh: int, nwl: float, ycut: float = 20.0):
        """Least-squares fit of u, v over |x - x0| <= nwl pi / alpha: coefficients, mask,
        relative rms residual of u and v (y <= ycut)."""
        m = np.abs(self.x - x0) <= nwl * np.pi / alpha
        X = self.x[m] - x0
        G = np.array([np.ones_like(X)] + [f(n * alpha * X) for n in range(1, nh + 1) for f in (np.cos, np.sin)]).T
        jy = self.y <= ycut
        coefs, err = {}, {}
        for name, q in (("u", self.u), ("v", self.v)):
            coefs[name] = np.linalg.lstsq(G, q[m], rcond=None)[0]
            r = q[m][:, jy] - G @ coefs[name][:, jy]
            err[name] = np.sqrt(np.mean(r ** 2) / np.mean(q[m][:, jy] ** 2))
        return coefs, err

    def alpha_at(self, x0: float, nh: int = 3, nwl: float = 2.0, method: str = "fit") -> float:
        """Local wavenumber. "hilbert": median over x0 +- 10 of the phase slope of v at the height
        of its maximum. "fit" (default): the alpha in 0.75..1.25 times that estimate that
        minimises the residual of the harmonic fit of u and v (the Hilbert phase slope oscillates
        when the wave has strong harmonics, and its short median can be off by 10 %)."""
        a0 = float(np.median(self.alpha_x[np.abs(self.x - x0) < 10]))
        if method == "hilbert":
            return a0
        J = lambda a: sum(self._fit(x0, a, nh, nwl)[1].values())
        grid = a0 * np.linspace(0.75, 1.25, 51)
        k = int(np.argmin([J(a) for a in grid]))
        lo, hi = grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]
        return float(minimize_scalar(J, bounds=(lo, hi), method="bounded", options=dict(xatol=1e-5)).x)

    def at(self, x0: float, nh: int = 3, nwl: float = 2.0, alpha_method: str = "fit") -> dict:
        """alpha, c, U_BF and harmonics u_n, v_n (n = -nh..nh) fitted over nwl wavelengths;
        fit_err: relative rms residual of the fit of u and v (y <= 20)."""
        alpha = self.alpha_at(x0, nh, nwl, alpha_method)
        coefs, err = self._fit(x0, alpha, nh, nwl)
        h = {}
        for name, coef in coefs.items():
            hn = {0: coef[0]}
            for n in range(1, nh + 1):
                hn[n] = 0.5 * (coef[2 * n - 1] - 1j * coef[2 * n])
                hn[-n] = np.conj(hn[n])
            h[name] = hn
        i0 = np.argmin(np.abs(self.x - x0))
        return dict(alpha=alpha, c=self.omega / alpha, u=h["u"], v=h["v"], ubf=self.ub[i0],
                    amp=2 * np.abs(h["u"][1]).max(), fit_err=err)

    def operator(self, x0: float, Re: float, n: int = 48, M: int = 3, ymax: float = 40.0,
                 yi: float = 2.5, nh: int = 3, scale: float = 1.0, power: bool = False,
                 alpha_method: str = "fit"):
        """SecondaryStability about the wave at x0, with the wave (not U_BF) scaled by `scale`.

        power=False: every harmonic and the mean-flow distortion u_0 are scaled by `scale`.
        power=True:  harmonic n is scaled by scale**|n| and u_0 by scale**2 (the weakly
                     nonlinear amplitude dependence), so a small wave tends to U_BF + scale u_1.
        """
        L = self.at(x0, nh, alpha_method=alpha_method)
        y, D1, D2 = grid(n, ymax, yi)
        spl = lambda p: _spline(self.y, p, y)
        f = (lambda k: scale ** (abs(k) if k else 2)) if power else (lambda k: scale)
        Un = {k: f(k) * spl(L["u"][k]) for k in L["u"] if k}
        Vn = {k: f(k) * spl(L["v"][k]) for k in L["v"] if k}
        Un[0] = spl(L["ubf"]) + f(0) * spl(L["u"][0]) - L["c"]
        Vn[0] = 0 * y
        return SecondaryStability(y, D1, D2, Un, Vn, L["alpha"], Re, M), L


def _spline(ysrc, p, y):
    f = lambda q: CubicSpline(ysrc, q)(np.minimum(y, ysrc[-1]))
    return f(p.real) + 1j * f(p.imag) if np.iscomplexobj(p) else f(p)
