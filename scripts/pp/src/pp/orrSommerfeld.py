"""
Local Orr-Sommerfeld analysis of a boundary layer profile.

This is a self contained Chebyshev collocation solver, written so that the
*eigenfunctions* are available -- the C++ solver in ``orrSommerfeldSolver``
only writes them for a single run, and the N factor synthesis in
``pp.localNfactor`` needs one solve per (station, frequency) pair, which is far
too many process launches.

Everything here works in the **global** non dimensionalisation of the DNS
(lengths in delta*_le, velocity in U_inf, Re = Re_delta*_le), so no rescaling
of the profile by the local delta* is needed and alpha, omega come out directly
in the units the N factor integral is taken in.

The two branches:

* ``temporal(alpha)`` -- real alpha, eigenvalue omega. Generalised problem of
  size n, cheap, and the input to a Gaster transformation.
* ``spatial(omega)``  -- real omega, eigenvalue alpha. The Orr-Sommerfeld
  operator is quartic in alpha, so this is a polynomial eigenvalue problem
  linearised onto a companion pencil of size 4n.

The spatial branch is what the N factor needs. Because a global 4n solve at
every (station, frequency) is expensive, ``spatial_newton`` refines a single
eigenvalue from a guess with a residual inverse iteration; the driver uses the
global solve only to seed the continuation and whenever the continuation fails.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.linalg as la
from numpy.polynomial import chebyshev as _cheb

__all__ = [
    "chebdif",
    "Mapping",
    "OrrSommerfeld",
    "gaster",
]


# --------------------------------------------------------------------------
# Chebyshev differentiation matrices
# --------------------------------------------------------------------------
def chebdif(n: int, m: int) -> tuple[np.ndarray, list[np.ndarray]]:
    """
    Chebyshev differentiation matrices of order 1..m on n Gauss-Lobatto nodes.

    Port of Weideman & Reddy's ``chebdif``. The higher order matrices are built
    by the recursion rather than by multiplying the first order one: powers of
    D1 amplify roundoff like n^(2k) and the fourth derivative of the
    Orr-Sommerfeld operator is exactly where that would bite.

    Returns ``(x, [D1, ..., Dm])`` with ``x`` descending from 1 to -1, which is
    the ordering the recursion is written for.
    """
    if n < 2:
        raise ValueError("need at least 2 nodes")

    ell = np.eye(n, dtype=bool)

    n1, n2 = n // 2, -(-n // 2)  # floor, ceil

    k = np.arange(n)
    th = k * np.pi / (n - 1)

    # nodes, written with sin() so that they stay symmetric to roundoff
    x = np.sin(np.pi * np.arange(n - 1, -n, -2) / (2 * (n - 1)))

    T = np.tile(th / 2, (n, 1)).T
    DX = 2.0 * np.sin(T.T + T) * np.sin(T.T - T)  # x[k] - x[j], by the identity
    DX = np.vstack((DX[:n1, :], -np.flipud(np.fliplr(DX[:n2, :]))))
    DX[ell] = 1.0

    C = la.toeplitz((-1.0) ** k)
    C[0, :] *= 2.0
    C[-1, :] *= 2.0
    C[:, 0] /= 2.0
    C[:, -1] /= 2.0

    Z = 1.0 / DX
    Z[ell] = 0.0

    D = np.eye(n)
    out: list[np.ndarray] = []
    for order in range(1, m + 1):
        D = order * Z * (C * np.tile(np.diag(D), (n, 1)).T - D)
        # negative sum trick: the diagonal is the one entry the recursion above
        # gets wrong, and row sums of a differentiation matrix must vanish
        D[ell] = -np.sum(D, axis=1)
        out.append(D.copy())

    return x, out


def clenshaw_curtis_weights(n: int) -> np.ndarray:
    """
    Clenshaw-Curtis weights for the n Gauss-Lobatto nodes of :func:`chebdif`,
    i.e. ``sum(w * f(x)) ~ int_{-1}^{1} f dx``, in the same descending order.
    """
    nn = n - 1
    c = np.zeros(n)
    k = np.arange(0, nn // 2 + 1)
    # standard FFT free construction
    theta = np.pi * np.arange(n) / nn
    w = np.zeros(n)
    for i in range(n):
        s = 0.0
        for j in k[1:]:
            s += np.cos(2 * j * theta[i]) / (4 * j * j - 1)
        w[i] = 1.0 - 2.0 * s
        if nn % 2 == 0:
            w[i] -= np.cos(nn * theta[i]) / (nn * nn - 1)
    w *= 2.0 / nn
    w[0] /= 2.0
    w[-1] /= 2.0
    del c
    return w


# --------------------------------------------------------------------------
# algebraic stretching of the wall normal coordinate
# --------------------------------------------------------------------------
@dataclass
class Mapping:
    """
    Malik's algebraic map from xi in [-1, 1] onto y in [ymin, ymax].

        y = ymin + a (1 + xi) / (b - xi)

    with ``a``, ``b`` fixed by requiring that half of the nodes sit below
    ``ymin + yi``. That clustering is what lets ~120 nodes resolve a
    Tollmien-Schlichting eigenfunction on a domain reaching far into the free
    stream.
    """

    ymin: float
    ymax: float
    yi: float

    def __post_init__(self) -> None:
        ly = self.ymax - self.ymin
        if not 0.0 < self.yi < 0.5 * ly:
            raise ValueError(
                f"need 0 < yi < (ymax - ymin)/2, got yi={self.yi}, L={ly}"
            )
        self.a = self.yi * ly / (ly - 2.0 * self.yi)
        self.b = 1.0 + 2.0 * self.a / ly

    def y(self, xi: np.ndarray) -> np.ndarray:
        return self.ymin + self.a * (1.0 + xi) / (self.b - xi)

    def dxi(self, y: np.ndarray) -> list[np.ndarray]:
        """d^k xi / dy^k for k = 1..4, evaluated at the nodes."""
        s = y - self.ymin + self.a
        f = self.a * (1.0 + self.b)
        return [
            f / s**2,
            -2.0 * f / s**3,
            6.0 * f / s**4,
            -24.0 * f / s**5,
        ]


def _chain_rule(
    dm: list[np.ndarray], a: list[np.ndarray]
) -> list[np.ndarray]:
    """
    Map D1..D4 in xi onto D1..D4 in y by Faa di Bruno, with
    ``a[k-1] = d^k xi / dy^k`` held as diagonals.

    Building D4 this way keeps the spectral accuracy of the xi matrices; taking
    the fourth power of the mapped first derivative matrix would not.
    """
    a1, a2, a3, a4 = (np.asarray(x) for x in a)
    d1, d2, d3, d4 = dm

    def dg(v: np.ndarray, m: np.ndarray) -> np.ndarray:
        return v[:, None] * m

    D1 = dg(a1, d1)
    D2 = dg(a1**2, d2) + dg(a2, d1)
    D3 = dg(a1**3, d3) + dg(3.0 * a1 * a2, d2) + dg(a3, d1)
    D4 = (
        dg(a1**4, d4)
        + dg(6.0 * a1**2 * a2, d3)
        + dg(4.0 * a1 * a3 + 3.0 * a2**2, d2)
        + dg(a4, d1)
    )
    return [D1, D2, D3, D4]


# --------------------------------------------------------------------------
# the solver
# --------------------------------------------------------------------------
class OrrSommerfeld:
    """
    Orr-Sommerfeld operator for one streamwise station.

    Parameters
    ----------
    y, u
        The base flow profile, as stored in the section files. ``y`` may start
        below zero (a station over the gap), in which case the no slip
        condition is applied on the gap floor.
    re
        Reynolds number in the *same* units as ``y`` -- the global Re of the
        DNS, not a local one, since no rescaling is done.
    n
        Number of collocation nodes.
    ymax
        Where to truncate the free stream. The eigenfunction has decayed like
        exp(-alpha y) long before the top of the section data, and pushing the
        truncation out costs resolution in the critical layer.
    yi
        Half of the nodes fall below ``ymin + yi``.
    """

    def __init__(
        self,
        y: np.ndarray,
        u: np.ndarray,
        re: float,
        n: int = 121,
        ymax: float = 45.0,
        yi: float | None = None,
        fit_degree: int = 40,
    ) -> None:
        y = np.asarray(y, dtype=float)
        u = np.asarray(u, dtype=float)
        order = np.argsort(y)
        y, u = y[order], u[order]
        # a profile sampled twice at the same y breaks the spline
        keep = np.concatenate(([True], np.diff(y) > 0.0))
        y, u = y[keep], u[keep]

        self.re = float(re)
        self.n = int(n)
        self.ymin = float(y[0])
        self.ymax = float(min(ymax, y[-1]))
        if yi is None:
            # clusters the nodes in the shear layer without starving the
            # critical layer, which for a TS mode sits at y = O(1)
            yi = min(4.0, 0.25 * (self.ymax - self.ymin))
        self.map = Mapping(self.ymin, self.ymax, float(yi))

        xi, dm = chebdif(self.n, 4)
        self.xi = xi
        self.yn = self.map.y(xi)
        self.D1, self.D2, self.D3, self.D4 = _chain_rule(dm, self.map.dxi(self.yn))

        # quadrature weights in y, for the L2 norm of the eigenfunction
        wxi = clenshaw_curtis_weights(self.n)
        dydxi = self.map.a * (1.0 + self.map.b) / (self.map.b - xi) ** 2
        self.w = wxi * dydxi

        # The profile is fitted by a truncated Chebyshev series in the mapped
        # coordinate rather than interpolated. The section files carry ~1e-4 of
        # noise in u (spectral element interpolation onto the sample points),
        # which is invisible in u itself but is ~1% of dU/dy and leaves U''
        # -- which enters the operator directly, next to the critical layer --
        # completely meaningless. Truncating the series at a degree where the
        # residual has reached that noise floor filters it out, and because the
        # degree stays below n the collocation D2 differentiates the fit
        # exactly, so U'' costs nothing extra and is consistent with the
        # discretisation of v''.
        inside = (y >= self.ymin) & (y <= self.ymax)
        yd, ud = y[inside], u[inside]
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
        coef, *_ = np.linalg.lstsq(V, ud, rcond=None)
        self.fit_rms = float(np.sqrt(np.mean((ud - V @ coef) ** 2)))
        self.U = _cheb.chebval(self.xi, coef)
        self.Upp = self.D2 @ self.U

        self.I = np.eye(self.n)
        # rows carrying the boundary conditions: v = 0 and v' = 0 at both ends.
        # Node 0 is y = ymax and node n-1 is y = ymin, because chebdif orders
        # xi descending.
        self._bc_rows = (0, 1, self.n - 2, self.n - 1)

    # -- boundary bordering ------------------------------------------------
    def _apply_bc(self, A: np.ndarray, B: np.ndarray) -> None:
        """Overwrite the four boundary rows of the pencil (A, B) in place."""
        r0, r1, rm2, rm1 = self._bc_rows
        for r in (r0, r1, rm2, rm1):
            A[r, :] = 0.0
            B[r, :] = 0.0
        A[r0, 0] = 1.0  # v(ymax) = 0
        A[rm1, -1] = 1.0  # v(ymin) = 0
        A[r1, :] = self.D1[0, :]  # v'(ymax) = 0
        A[rm2, :] = self.D1[-1, :]  # v'(ymin) = 0

    def _bc_block(self) -> np.ndarray:
        """The 4 x n matrix of boundary condition rows."""
        out = np.zeros((4, self.n))
        out[0, 0] = 1.0
        out[1, :] = self.D1[0, :]
        out[2, :] = self.D1[-1, :]
        out[3, -1] = 1.0
        return out

    # -- temporal branch ---------------------------------------------------
    def temporal(self, alpha: float) -> tuple[np.ndarray, np.ndarray]:
        """
        Eigenvalues ``omega`` and eigenfunctions ``v`` at fixed real ``alpha``.

            omega (D2 - a^2) v = [ a (U (D2 - a^2) - U'') + i/Re (D2 - a^2)^2 ] v
        """
        a2 = alpha * alpha
        L = self.D2 - a2 * self.I
        A = (
            alpha * (self.U[:, None] * L - np.diag(self.Upp))
            + 1j / self.re * (L @ L)
        )
        B = L.astype(complex)
        A = A.astype(complex)
        self._apply_bc(A, B)

        w, v = la.eig(A, B)
        finite = np.isfinite(w) & (np.abs(w) < 1e6)
        return w[finite], v[:, finite]

    # -- spatial branch ----------------------------------------------------
    def _spatial_coeffs(self, omega: complex) -> list[np.ndarray]:
        """
        Coefficients Q0, Q1, Q2 of the quadratic ``sum_k Q_k alpha^k z = 0``
        in the split variables ``z = (v, q)``, ``q = (D2 - alpha^2) v``.

        Collecting powers of alpha in the Orr-Sommerfeld equation directly
        gives a *quartic* in alpha, whose companion pencil mixes D4 with the
        identity -- entries spanning n^8, which destroys the small physical
        root. Splitting off q (Bridges & Morris 1984) drops the highest
        derivative to D2 and makes the same sized pencil well behaved.

            (i a U - i w) q - i a U'' v = 1/Re (D2 - a^2) q
            q = (D2 - a^2) v

        The four boundary rows are written into Q0 with Q1 = Q2 = 0 there, so
        that v = v' = 0 holds independently of alpha.
        """
        n = self.n
        Z = np.zeros((n, n), dtype=complex)
        I = self.I.astype(complex)

        Q2 = np.block([[Z, I / self.re], [-I, Z]])
        Q1 = np.block(
            [[-1j * np.diag(self.Upp), 1j * np.diag(self.U)], [Z, Z]]
        )
        Q0 = np.block(
            [
                [Z, -1j * omega * I - self.D2.astype(complex) / self.re],
                [self.D2.astype(complex), -I],
            ]
        )

        # boundary bordering: v' = 0 replaces the two end rows of the momentum
        # equation, v = 0 the two end rows of the definition of q
        bc = self._bc_block()
        rows = (0, n - 1, n + 0, n + n - 1)
        for M in (Q0, Q1, Q2):
            for r in rows:
                M[r, :] = 0.0
        Q0[0, :n] = bc[1]  # v'(ymax) = 0
        Q0[n - 1, :n] = bc[2]  # v'(ymin) = 0
        Q0[n + 0, :n] = bc[0]  # v(ymax)  = 0
        Q0[n + n - 1, :n] = bc[3]  # v(ymin)  = 0

        return [Q0, Q1, Q2]

    def spatial(
        self, omega: complex, vectors: bool = True
    ) -> tuple[np.ndarray, np.ndarray | None]:
        """
        Eigenvalues ``alpha`` and eigenfunctions ``v`` at fixed ``omega``.

        Linearised on the companion pencil ``A y = alpha B y`` with
        ``y = (z, alpha z)``. The four boundary rows make Q2 singular, so the
        infinite eigenvalues they produce are filtered on the way out.

        ``vectors=False`` skips the eigenvectors, which roughly halves the
        cost. The N factor never needs them -- normalising every mode to unit
        L2 norm is exactly what removes the eigenfunction from the answer --
        and the mode is selected on alpha alone.
        """
        n = self.n
        Q0, Q1, Q2 = self._spatial_coeffs(omega)
        m = 2 * n
        I2 = np.eye(m, dtype=complex)
        Zm = np.zeros((m, m), dtype=complex)

        A = np.block([[Zm, I2], [-Q0, -Q1]])
        B = np.block([[I2, Zm], [Zm, Q2]])

        if not vectors:
            w = la.eig(A, B, right=False)
            good = np.isfinite(w) & (np.abs(w) < 1e4)
            return w[good], None

        w, v = la.eig(A, B)
        good = np.isfinite(w) & (np.abs(w) < 1e4)
        return w[good], v[:n, good]

    def spatial_newton(
        self,
        omega: complex,
        alpha0: complex,
        v0: np.ndarray | None = None,
        tol: float = 1e-11,
        maxiter: int = 40,
        trust: float | None = None,
    ) -> tuple[complex, np.ndarray] | None:
        """
        Refine one spatial eigenvalue from a guess by residual inverse
        iteration on the quadratic.

        Each step is a 2n x 2n solve instead of a 4n x 4n eigendecomposition,
        which is what makes a frequency by station sweep affordable.

        ``v0`` is the eigenfunction the guess belongs to, from the neighbouring
        station or frequency. It matters: inverse iteration converges to
        whatever the starting vector has most of, so without it the iteration
        happily lands on a different branch of the spectrum even when
        ``alpha0`` is good. ``trust`` rejects a result that has wandered
        further than that from the guess.

        Returns ``None`` when the iteration does not settle or leaves the trust
        radius, and the caller then falls back to :meth:`spatial`.
        """
        Q0, Q1, Q2 = self._spatial_coeffs(omega)

        def P(a: complex) -> np.ndarray:
            return Q0 + a * Q1 + a * a * Q2

        def dP(a: complex) -> np.ndarray:
            return Q1 + 2.0 * a * Q2

        a = complex(alpha0)
        n = self.n
        z = np.zeros(2 * n, dtype=complex)
        if v0 is not None:
            v0 = np.asarray(v0, dtype=complex)
            z[:n] = v0
            z[n:] = (self.D2 - a * a * self.I) @ v0
        else:
            decay = np.exp(-(self.yn - self.ymin) / max(0.05 * self.ymax, 1e-3))
            z[:n] = decay
            z[n:] = decay
        nz = la.norm(z)
        if nz == 0.0 or not np.isfinite(nz):
            return None
        z /= nz

        for _ in range(maxiter):
            try:
                lu = la.lu_factor(P(a))
                w = la.lu_solve(lu, dP(a) @ z)
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
                return a, z[:n]
        return None

    # -- helpers -----------------------------------------------------------
    def velocity(self, v: np.ndarray, alpha: complex) -> tuple[np.ndarray, np.ndarray]:
        """(u_hat, v_hat) from the wall normal eigenfunction, via continuity."""
        return 1j * (self.D1 @ v) / alpha, v

    def norm(self, v: np.ndarray, alpha: complex) -> float:
        """sqrt(int (|u|^2 + |v|^2) dy), the quantity the DNS L2 amplitude is."""
        uh, vh = self.velocity(v, alpha)
        return float(
            np.sqrt(np.abs(np.sum(self.w * (np.abs(uh) ** 2 + np.abs(vh) ** 2))))
        )


# --------------------------------------------------------------------------
# Gaster transformation (kept here so the spatial branch can be cross checked)
# --------------------------------------------------------------------------
def gaster(
    omega_r: np.ndarray, omega_i: np.ndarray, alpha_r: np.ndarray, order: int = 2
) -> tuple[np.ndarray, np.ndarray]:
    """
    Map a temporal sweep omega(alpha_r) onto the spatial branch.

    Same expansion as ``scripts/getlocalGrowthRateOS.py``: order 1 is
    alpha_i = -omega_i / c_g, order 2 keeps the curvature of the temporal
    dispersion relation. Returns ``(alpha_i, omega)``.
    """
    cg = np.abs(np.gradient(omega_r, alpha_r, edge_order=2))
    with np.errstate(divide="ignore", invalid="ignore"):
        ai1 = -omega_i / cg
    if order == 1:
        return ai1, omega_r.copy()

    h1 = alpha_r[1:-1] - alpha_r[:-2]
    h2 = alpha_r[2:] - alpha_r[1:-1]

    def d2(f: np.ndarray) -> np.ndarray:
        inner = (
            2.0
            * (h2 * f[:-2] - (h1 + h2) * f[1:-1] + h1 * f[2:])
            / (h1 * h2 * (h1 + h2))
        )
        return np.concatenate(([inner[0]], inner, [inner[-1]]))

    doi = np.gradient(omega_i, alpha_r, edge_order=2)
    d2oi, d2or = d2(omega_i), d2(omega_r)

    disc = cg**2 + 2.0 * d2oi * omega_i
    denom = cg + np.sqrt(np.where(disc > 0.0, disc, 0.0))
    usable = (disc > 0.0) & (denom != 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ai = np.where(usable, -2.0 * omega_i / denom, ai1)
    om = np.where(usable, omega_r - ai * doi - 0.5 * ai**2 * d2or, omega_r)
    return ai, om
