"""
Local (parallel flow) compressible linear stability of a boundary layer.

The compressible counterpart of :mod:`pp.orrSommerfeld`: same Chebyshev
collocation, same Malik stretching, same profile fitting, same temporal /
spatial / Newton interface, so it can be dropped into the N factor machinery.

Equations
---------
Compressible Navier-Stokes for an ideal gas, Stokes hypothesis
(lambda = -2/3 mu), constant Prandtl number, viscosity from Sutherland's law
(as in Nektar++ with ``ViscosityType = Variable``), non dimensionalised with the
free stream values rho_inf, U_inf, T_inf, mu_inf and the length of the profile
(delta*_le for the DNS sections), so that

    p = rho T / (gamma Ma^2),    Re = rho_inf U_inf L / mu_inf.

The base flow is parallel, U(y) and T(y) with V = W = 0, and the pressure is
uniform, so rho = 1 / T. Perturbations are normal modes

    q'(x, y, z, t) = q(y) exp(i (alpha x + beta z - omega t))

in the primitive variables ``q = (u, v, w, T, p)``. Pressure, not density, is
the thermodynamic unknown (as in Malik 1990): the continuity equation is then
regular as Ma -> 0, and at Ma = 0 with T = 1 the system reduces *exactly* to
Orr-Sommerfeld plus a decoupled temperature equation, which is one of the
checks in ``tests/testCompressibleStability.py``.

The linear operator is

    L(alpha, omega) q = (L0 + alpha L1 + alpha^2 L2 + omega Lw) q = 0,

i.e. at most quadratic in alpha (second derivatives only, no fourth order
operator to split) and linear in omega. Hence

* ``temporal(alpha)`` -- generalised eigenproblem of size 5n,
* ``spatial(omega)``  -- quadratic eigenproblem, companion pencil of size 10n,
* ``spatial_newton``  -- residual inverse iteration on the quadratic, size 5n.

Boundary conditions: u = v = w = 0 at the wall and in the free stream;
T = 0 at both ends (``thermal_bc="isothermal"``, the usual choice, and the
right one for the Nektar++ ``WallViscous`` + ``Twall`` set up) or dT/dy = 0 at
the wall (``"adiabatic"``). The pressure carries no boundary condition: the
continuity equation is collocated at every node, including the boundaries.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import scipy.linalg as la
from numpy.polynomial import chebyshev as _cheb
from scipy.integrate import cumulative_simpson, simpson, solve_bvp

from pp.orrSommerfeld import Mapping, _chain_rule, chebdif, clenshaw_curtis_weights

__all__ = [
    "Gas",
    "CompressibleStability",
    "compressible_blasius",
    "displacement_thickness",
    "VARS",
]

#: ordering of the unknowns inside one eigenvector
VARS = ("u", "v", "w", "T", "p")
_U, _V, _W, _T, _P = range(5)


# --------------------------------------------------------------------------
# thermodynamics
# --------------------------------------------------------------------------
@dataclass
class Gas:
    """
    Fluid properties, all dimensionless except ``tinf`` and ``sutherland``.

    The defaults are those of the Nektar++ session files: gamma = 1.4,
    Pr = 0.72, Sutherland with S = 110 K and T_inf = 225 K (``Tref`` in the
    session). ``viscosity="constant"`` gives mu = 1 everywhere, which is useful
    for comparisons with the incompressible solver.
    """

    mach: float
    gamma: float = 1.4
    pr: float = 0.72
    tinf: float = 225.0
    sutherland: float = 110.0
    viscosity: str | Callable[[np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray]] = "sutherland"

    def mu(self, T: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """mu, dmu/dT, d2mu/dT2 at the (dimensionless) temperature T."""
        T = np.asarray(T, dtype=float)
        if callable(self.viscosity):
            return self.viscosity(T)
        if self.viscosity == "constant":
            one = np.ones_like(T)
            return one, 0.0 * one, 0.0 * one
        if self.viscosity != "sutherland":
            raise ValueError(f"unknown viscosity law {self.viscosity!r}")
        s = self.sutherland / self.tinf
        mu = T**1.5 * (1.0 + s) / (T + s)
        g = 1.5 / T - 1.0 / (T + s)  # dlog(mu)/dT
        mu1 = mu * g
        mu2 = mu * (g * g - 1.5 / T**2 + 1.0 / (T + s) ** 2)
        return mu, mu1, mu2

    def temperature(self, rho, rhou, E, rhov=None) -> np.ndarray:
        """T from the conservative variables of the Nektar++ fields."""
        rho = np.asarray(rho, dtype=float)
        ke = 0.5 * np.asarray(rhou) ** 2 / rho
        if rhov is not None:
            ke = ke + 0.5 * np.asarray(rhov) ** 2 / rho
        p = (self.gamma - 1.0) * (np.asarray(E) - ke)
        return self.gamma * self.mach**2 * p / rho


# --------------------------------------------------------------------------
# compressible similarity solution (for tests and for quick studies)
# --------------------------------------------------------------------------
def compressible_blasius(
    gas: Gas,
    twall: float | None = 1.0,
    eta_max: float = 20.0,
    npts: int = 8001,
    y_far: float = 400.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Self similar zero pressure gradient compressible boundary layer.

    Same equations as ``blasius/solve_bvp.py`` (Howarth-Dorodnitsyn variable,
    Chapman-Rubesin parameter C = rho mu = mu(T)/T):

        (C f'')' + f f'' / 2 = 0
        (C g' / Pr)' + f g' / 2 + (gamma - 1) Ma^2 C f''^2 = 0,   g = T

    ``twall=None`` gives an adiabatic wall, a number an isothermal one
    (``twall=1`` is the ``Twall = TInf`` of the session files).

    The similarity problem is solved up to ``eta_max`` and the profile is then
    padded with free stream values up to ``y = y_far``: the solver silently
    truncates its domain at the end of the profile, and a TS eigenfunction at
    small alpha needs a far field boundary many boundary layer thicknesses out.

    Returns ``(y, U, T)`` with y in Blasius lengths ``sqrt(nu_inf x / U_inf)``,
    i.e. ``Re = sqrt(Re_x)`` in the stability solver. Divide y by
    :func:`displacement_thickness` to work in delta* units instead.
    """

    def C(g):
        mu, mu1, _ = gas.mu(g)
        return mu / g, (mu1 / g - mu / g**2)

    def odes(eta, Y):
        f, f1, f2, g, g1 = Y
        c, cT = C(g)
        dc = cT * g1
        df2 = -f2 * (dc + 0.5 * f) / c
        dg1 = (
            -g1 * (dc + 0.5 * gas.pr * f) / c
            - (gas.gamma - 1.0) * gas.pr * gas.mach**2 * f2**2
        )
        return np.vstack((f1, f2, df2, g1, dg1))

    def bc(ya, yb):
        wall = ya[4] if twall is None else ya[3] - twall
        return np.array([ya[0], ya[1], yb[1] - 1.0, yb[3] - 1.0, wall])

    eta = np.linspace(0.0, eta_max, 401)
    taw = 1.0 + np.sqrt(gas.pr) * 0.5 * (gas.gamma - 1.0) * gas.mach**2
    tw0 = taw if twall is None else twall
    guess = np.zeros((5, eta.size))
    guess[1] = np.tanh(0.5 * eta)
    guess[0] = 2.0 * np.log(np.cosh(0.5 * eta))
    guess[2] = 0.5 / np.cosh(0.5 * eta) ** 2
    guess[3] = 1.0 + (tw0 - 1.0) * np.exp(-0.5 * eta**2 / 4.0)
    guess[4] = -(tw0 - 1.0) * eta / 4.0 * np.exp(-0.5 * eta**2 / 4.0)
    sol = solve_bvp(odes, bc, eta, guess, tol=1e-10, max_nodes=200000)
    if not sol.success:
        raise RuntimeError(f"similarity solution failed: {sol.message}")

    eta = np.linspace(0.0, eta_max, npts)
    f, f1, f2, g, g1 = sol.sol(eta)
    # physical normal coordinate: dy = (rho_inf / rho) d eta = T d eta
    y = cumulative_simpson(g, x=eta, initial=0.0)
    if y_far > y[-1]:
        pad = np.linspace(y[-1], y_far, int(np.ceil((y_far - y[-1]) / 0.1)) + 1)[1:]
        y = np.concatenate((y, pad))
        f1 = np.concatenate((f1, np.ones_like(pad)))
        g = np.concatenate((g, np.ones_like(pad)))
    return y, f1, g


def displacement_thickness(y: np.ndarray, u: np.ndarray, T: np.ndarray) -> float:
    """delta* = int (1 - rho u) dy, with rho = 1/T."""
    return float(simpson(1.0 - u / T, x=y))


# --------------------------------------------------------------------------
# the solver
# --------------------------------------------------------------------------
class CompressibleStability:
    """
    Compressible linear stability operator for one streamwise station.

    Parameters
    ----------
    y, u, T
        Base flow profile. ``T`` is the dimensionless temperature
        T / T_inf; for DNS sections use :meth:`Gas.temperature` (or
        :meth:`from_conservative`). ``y`` may start below zero (station over a
        gap): the wall conditions go on the lowest point.
    re
        Reynolds number in the units of ``y`` (Re_delta*_le for the DNS).
    gas
        :class:`Gas` (Mach number, gamma, Pr, viscosity law).
    beta
        Spanwise wavenumber (0 for 2D waves).
    n, ymax, yi, fit_degree
        As in :class:`pp.orrSommerfeld.OrrSommerfeld`.
    thermal_bc
        ``"isothermal"`` (T' = 0 at the wall, default) or ``"adiabatic"``
        (dT'/dy = 0). In the free stream T' = 0 always.
    """

    def __init__(
        self,
        y: np.ndarray,
        u: np.ndarray,
        T: np.ndarray,
        re: float,
        gas: Gas,
        beta: float = 0.0,
        n: int = 121,
        ymax: float = 45.0,
        yi: float | None = None,
        fit_degree: int = 40,
        thermal_bc: str = "isothermal",
    ) -> None:
        y = np.asarray(y, dtype=float)
        u = np.asarray(u, dtype=float)
        T = np.asarray(T, dtype=float)
        order = np.argsort(y)
        y, u, T = y[order], u[order], T[order]
        keep = np.concatenate(([True], np.diff(y) > 0.0))
        y, u, T = y[keep], u[keep], T[keep]

        if thermal_bc not in ("isothermal", "adiabatic"):
            raise ValueError("thermal_bc must be 'isothermal' or 'adiabatic'")
        self.thermal_bc = thermal_bc
        self.re = float(re)
        self.gas = gas
        self.beta = float(beta)
        self.n = int(n)
        self.ymin = float(y[0])
        self.ymax = float(min(ymax, y[-1]))
        if yi is None:
            yi = min(4.0, 0.25 * (self.ymax - self.ymin))
        self.map = Mapping(self.ymin, self.ymax, float(yi))

        xi, dm = chebdif(self.n, 2)
        self.xi = xi
        self.yn = self.map.y(xi)
        # only D1, D2 are needed; pad the xi list so _chain_rule can be reused
        z = np.zeros_like(dm[0])
        self.D1, self.D2, _, _ = _chain_rule(dm + [z, z], self.map.dxi(self.yn))

        wxi = clenshaw_curtis_weights(self.n)
        dydxi = self.map.a * (1.0 + self.map.b) / (self.map.b - xi) ** 2
        self.w = wxi * dydxi

        # base flow by a truncated Chebyshev fit in the mapped coordinate: see
        # OrrSommerfeld.__init__ for why (noise in the section files)
        inside = (y >= self.ymin) & (y <= self.ymax)
        yd = y[inside]
        if fit_degree >= self.n - 1:
            raise ValueError("fit_degree must stay below n - 1 for D2 to be exact")
        if yd.size <= fit_degree:
            raise ValueError(
                f"profile has {yd.size} points below ymax = {self.ymax}, "
                f"too few for a degree {fit_degree} fit"
            )
        sd = yd - self.ymin
        xid = (sd * self.map.b - self.map.a) / (sd + self.map.a)
        V = _cheb.chebvander(xid, fit_degree)
        self.fit_rms = {}
        fitted = {}
        for name, f in (("u", u[inside]), ("T", T[inside])):
            coef, *_ = np.linalg.lstsq(V, f, rcond=None)
            self.fit_rms[name] = float(np.sqrt(np.mean((f - V @ coef) ** 2)))
            fitted[name] = _cheb.chebval(self.xi, coef)
        self.U, self.Tb = fitted["u"], fitted["T"]

        # A least squares fit is only pinned down where there is data: if the
        # profile is sampled more coarsely than the clustered collocation nodes
        # (near the wall, typically), a high degree fit can oscillate between
        # samples while its residual at the samples stays tiny. Compare with
        # plain interpolation at the nodes to catch that.
        self.fit_dev = {}
        for name, f, fit in (("u", u, self.U), ("T", T, self.Tb)):
            ref = np.interp(self.yn, y, f)
            self.fit_dev[name] = float(np.max(np.abs(fit - ref)))
        if np.any(self.Tb <= 0.0) or max(self.fit_dev.values()) > 1e-2:
            raise ValueError(
                "the Chebyshev fit of the base flow departs from the data "
                f"between samples (max deviation {self.fit_dev}); lower "
                "fit_degree or sample the profile more finely"
            )
        self.Up, self.Upp = self.D1 @ self.U, self.D2 @ self.U
        self.Tp, self.Tpp = self.D1 @ self.Tb, self.D2 @ self.Tb

        self.I = np.eye(self.n)
        self._build()

    @classmethod
    def from_conservative(
        cls, y, rho, rhou, E, re, gas: Gas, rhov=None, **kw
    ) -> "CompressibleStability":
        """Build from the Nektar++ conservative fields (rho, rhou, [rhov], E)."""
        rho = np.asarray(rho, dtype=float)
        u = np.asarray(rhou, dtype=float) / rho
        T = gas.temperature(rho, rhou, E, rhov)
        return cls(y, u, T, re, gas, **kw)

    # -- operator ------------------------------------------------------------
    def _build(self) -> None:
        """
        Assemble L0, L1, L2, Lw (each 5n x 5n) with boundary rows already
        bordered. Row block = equation, column block = variable, both in the
        order of ``VARS``; the equations are x, y, z momentum, energy,
        continuity.
        """
        n, g = self.n, self.gas
        gam, M2, Pr, Re, b = g.gamma, g.mach**2, g.pr, self.re, self.beta
        D1, D2, I = self.D1, self.D2, self.I

        T, Tp, Tpp = self.Tb, self.Tp, self.Tpp
        U, Up, Upp = self.U, self.Up, self.Upp
        rho = 1.0 / T
        mu, muT, muTT = g.mu(T)
        lam, lamT = -2.0 / 3.0 * mu, -2.0 / 3.0 * muT
        mup, lamp = muT * Tp, lamT * Tp  # d/dy of the base mu, lambda
        dmuUp = muTT * Tp * Up + muT * Upp  # d/dy (mu_T U')
        dmuTp = muTT * Tp**2 + muT * Tpp  # d/dy (mu_T T')

        def d(f):
            return np.diag(np.asarray(f, dtype=complex) * np.ones(n))

        Z = np.zeros((n, n))
        L0 = [[Z.astype(complex) for _ in range(5)] for _ in range(5)]
        L1 = [[Z.astype(complex) for _ in range(5)] for _ in range(5)]
        L2 = [[Z.astype(complex) for _ in range(5)] for _ in range(5)]
        Lw = [[Z.astype(complex) for _ in range(5)] for _ in range(5)]
        u, v, w, t, p = _U, _V, _W, _T, _P
        ec = 5 - 1  # continuity row block
        ee = _T  # energy row block
        k = gam * (gam - 1.0) * M2 / Re  # viscous heating
        kc = gam / (Re * Pr)  # conduction

        # x momentum ---------------------------------------------------------
        L0[u][u] = -(d(mup) @ D1 + d(mu) @ D2) / Re + d(b * b * mu / Re)
        L1[u][u] = d(1j * rho * U)
        L2[u][u] = d((2.0 * mu + lam) / Re)
        Lw[u][u] = d(-1j * rho)
        L0[u][v] = d(rho * Up)
        L1[u][v] = -1j / Re * (d(lam + mu) @ D1 + d(mup))
        L1[u][w] = d(b * (lam + mu) / Re)
        L0[u][t] = -(d(dmuUp) + d(muT * Up) @ D1) / Re
        L1[u][p] = 1j * I

        # y momentum ---------------------------------------------------------
        L1[v][u] = -1j / Re * (d(mu + lam) @ D1 + d(lamp))
        L0[v][v] = (
            -(d(2.0 * mu + lam) @ D2 + d(2.0 * mup + lamp) @ D1) / Re
            + d(b * b * mu / Re)
        )
        L1[v][v] = d(1j * rho * U)
        L2[v][v] = d(mu / Re)
        Lw[v][v] = d(-1j * rho)
        L0[v][w] = -1j * b / Re * (d(lam + mu) @ D1 + d(lamp))
        L1[v][t] = d(-1j * muT * Up / Re)
        L0[v][p] = D1.astype(complex)

        # z momentum ---------------------------------------------------------
        L1[w][u] = d(b * (mu + lam) / Re)
        L0[w][v] = -1j * b / Re * (d(mu + lam) @ D1 + d(mup))
        L0[w][w] = -(d(mup) @ D1 + d(mu) @ D2) / Re + d(b * b * (2.0 * mu + lam) / Re)
        L1[w][w] = d(1j * rho * U)
        L2[w][w] = d(mu / Re)
        Lw[w][w] = d(-1j * rho)
        L0[w][p] = 1j * b * I

        # energy ---------------------------------------------------------------
        L0[ee][u] = -2.0 * k * d(mu * Up) @ D1
        L1[ee][u] = 1j * (gam - 1.0) * I
        L0[ee][v] = d(rho * Tp) + (gam - 1.0) * D1
        L1[ee][v] = d(-2j * k * mu * Up)
        L0[ee][w] = 1j * (gam - 1.0) * b * I
        L0[ee][t] = (
            -kc * (d(mu) @ D2 + d(2.0 * muT * Tp) @ D1 + d(dmuTp - b * b * mu))
            - d(k * muT * Up**2)
        )
        L1[ee][t] = d(1j * rho * U)
        L2[ee][t] = d(kc * mu)
        Lw[ee][t] = d(-1j * rho)

        # continuity (divided by rho) -----------------------------------------
        L1[ec][u] = 1j * I
        L0[ec][v] = D1 - d(Tp / T)
        L0[ec][w] = 1j * b * I
        L1[ec][t] = d(-1j * U / T)
        Lw[ec][t] = d(1j / T)
        L1[ec][p] = d(1j * gam * M2 * U)
        Lw[ec][p] = -1j * gam * M2 * I

        mats = [np.block(X) for X in (L0, L1, L2, Lw)]

        # boundary bordering: node 0 is ymax, node n-1 the wall ---------------
        rows: list[tuple[int, np.ndarray]] = []
        for var in (u, v, w):
            for node in (0, n - 1):
                e = np.zeros(5 * n)
                e[var * n + node] = 1.0
                rows.append((var * n + node, e))
        e = np.zeros(5 * n)
        e[t * n + 0] = 1.0
        rows.append((t * n + 0, e))
        e = np.zeros(5 * n)
        if self.thermal_bc == "isothermal":
            e[t * n + n - 1] = 1.0
        else:
            e[t * n : (t + 1) * n] = D1[-1, :]
        rows.append((t * n + n - 1, e))
        for M in mats:
            for r, _ in rows:
                M[r, :] = 0.0
        for r, e in rows:
            mats[0][r, :] = e

        self.L0, self.L1, self.L2, self.Lw = mats

        # For 2D waves w decouples (its equation and every coupling to it are
        # proportional to beta), so it is dropped from the eigenproblems.
        keep = [x for x in range(5) if not (x == _W and self.beta == 0.0)]
        self._idx = np.concatenate([np.arange(x * n, (x + 1) * n) for x in keep])
        ix = np.ix_(self._idx, self._idx)
        self._L0, self._L1, self._L2, self._Lw = (M[ix] for M in mats)
        # unknowns that actually carry alpha^2: only these need an alpha q copy
        # in the companion pencil (the pressure, and the boundary values, which
        # are fixed by the bordering, do not)
        self._sq = np.flatnonzero(np.any(self._L2 != 0.0, axis=0))

    def _expand(self, q: np.ndarray) -> np.ndarray:
        """Reduced eigenvector(s) -> full 5n layout of ``VARS``."""
        out = np.zeros((5 * self.n,) + q.shape[1:], dtype=complex)
        out[self._idx] = q
        return out

    def operator(self, alpha: complex, omega: complex) -> np.ndarray:
        """L(alpha, omega) on the full 5n unknowns; its null vector is the mode."""
        return self.L0 + alpha * self.L1 + alpha * alpha * self.L2 + omega * self.Lw

    # -- temporal branch ---------------------------------------------------
    def temporal(self, alpha: complex, vectors: bool = True):
        """Eigenvalues ``omega`` (and 5n eigenvectors) at fixed ``alpha``."""
        A = self._L0 + alpha * self._L1 + alpha * alpha * self._L2
        B = -self._Lw
        if not vectors:
            w = la.eig(A, B, right=False)
            good = np.isfinite(w) & (np.abs(w) < 1e6)
            return w[good], None
        w, q = la.eig(A, B)
        good = np.isfinite(w) & (np.abs(w) < 1e6)
        return w[good], self._expand(q[:, good])

    # -- spatial branch ----------------------------------------------------
    def _spatial_coeffs(self, omega: complex) -> list[np.ndarray]:
        return [self._L0 + omega * self._Lw, self._L1, self._L2]

    def spatial(self, omega: complex, vectors: bool = True):
        """
        Eigenvalues ``alpha`` (and 5n eigenvectors) at fixed ``omega``.

        Linearisation of ``(Q0 + alpha Q1 + alpha^2 Q2) q = 0`` with
        ``r = alpha q_S`` kept only for the unknowns S that Q2 acts on:

            [ -Q0  0 ] [q]         [ Q1   Q2[:, S] ] [q]
            [  0   I ] [r] = alpha [ E_S     0     ] [r]

        For 2D waves this is a 7n pencil instead of the 10n of the plain
        companion form, about 3x cheaper. The rows without alpha^2 terms still
        produce infinite eigenvalues, which are filtered.
        """
        Q0, Q1, Q2 = self._spatial_coeffs(omega)
        m, S = Q0.shape[0], self._sq
        k = S.size
        A = np.zeros((m + k, m + k), dtype=complex)
        B = np.zeros((m + k, m + k), dtype=complex)
        A[:m, :m] = -Q0
        A[m:, m:] = np.eye(k)
        B[:m, :m] = Q1
        B[:m, m:] = Q2[:, S]
        B[m + np.arange(k), S] = 1.0
        if not vectors:
            w = la.eig(A, B, right=False)
            good = np.isfinite(w) & (np.abs(w) < 1e4)
            return w[good], None
        w, z = la.eig(A, B)
        good = np.isfinite(w) & (np.abs(w) < 1e4)
        return w[good], self._expand(z[:m, good])

    def spatial_newton(
        self,
        omega: complex,
        alpha0: complex,
        q0: np.ndarray | None = None,
        tol: float = 1e-11,
        maxiter: int = 40,
        trust: float | None = None,
    ) -> tuple[complex, np.ndarray] | None:
        """
        Refine one spatial eigenvalue from a guess, by residual inverse
        iteration on the quadratic (see ``OrrSommerfeld.spatial_newton``).
        ``q0`` is the full 5n eigenvector of the neighbouring solution; the
        returned eigenvector is 5n as well.
        """
        Q0, Q1, Q2 = self._spatial_coeffs(omega)
        a = complex(alpha0)
        if q0 is not None:
            z = np.asarray(q0, dtype=complex)[self._idx].copy()
        else:
            decay = np.exp(-(self.yn - self.ymin) / max(0.05 * self.ymax, 1e-3))
            z = np.tile(decay, self._idx.size // self.n).astype(complex)
        nz = la.norm(z)
        if nz == 0.0 or not np.isfinite(nz):
            return None
        z /= nz
        for _ in range(maxiter):
            try:
                lu = la.lu_factor(Q0 + a * Q1 + a * a * Q2)
                w = la.lu_solve(lu, (Q1 + 2.0 * a * Q2) @ z)
            except (la.LinAlgError, ValueError):
                return None
            denom = np.vdot(z, w)
            if denom == 0.0 or not np.isfinite(denom):
                return None
            a_new = a - 1.0 / denom
            nw = la.norm(w)
            if nw == 0.0 or not np.isfinite(nw) or not np.isfinite(a_new):
                return None
            z = w / nw
            converged = abs(a_new - a) < tol * max(1.0, abs(a_new))
            a = a_new
            if converged:
                if trust is not None and abs(a - alpha0) > trust:
                    return None
                return a, self._expand(z)
        return None

    # -- helpers -----------------------------------------------------------
    def split(self, q: np.ndarray) -> dict[str, np.ndarray]:
        """Eigenvector -> {'u', 'v', 'w', 'T', 'p', 'rho'} profiles."""
        n = self.n
        out = {name: q[i * n : (i + 1) * n] for i, name in enumerate(VARS)}
        g = self.gas
        rho = 1.0 / self.Tb
        out["rho"] = rho * (g.gamma * g.mach**2 * out["p"] - rho * out["T"])
        return out

    def norm(self, q: np.ndarray, components=("u", "v")) -> float:
        """sqrt(int sum |q_c|^2 dy) over the chosen components."""
        f = self.split(q)
        s = sum(np.abs(f[c]) ** 2 for c in components)
        return float(np.sqrt(np.abs(np.sum(self.w * s))))
