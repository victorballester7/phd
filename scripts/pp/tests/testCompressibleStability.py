"""
Benchmarks for the compressible local stability solver in
``pp.compressibleStability``.

Run directly:  python tests/testCompressibleStability.py

Each check compares against a number that does not depend on this code:

1. The operator against a symbolic linearisation (sympy) of the full
   compressible Navier-Stokes equations, at Ma = 0.8, Sutherland viscosity, a
   non uniform temperature profile, an oblique wave and complex alpha, omega.
   This checks every coefficient of the operator, independently of the
   eigenvalue solvers.
2. Ma = 0, T = 1: the system must reduce exactly to Orr-Sommerfeld, so it is
   compared with ``pp.orrSommerfeld`` (temporal and spatial).
3. Incompressible Blasius, spatial: Jordinson (1970), Re_delta* = 998,
   omega = 0.1122 -> alpha = 0.3086 - 0.0057i (as tabulated in Danabasoglu &
   Biringen, NASA CR-4247, Table 2, whose own finite difference value
   0.30870 - 0.00564i is less accurate: the Orr-Sommerfeld solver gives the
   same as this code, 0.308591 - 0.005708i).
4. Malik (1990), J. Comput. Phys. 86, case 1: Ma = 0.5, adiabatic plate,
   T0 = 278 K, R = 2000, alpha = 0.1 -> omega = 0.0290817 + 0.0022441i
   (temporal, first mode).
5. Malik (1990) case 3: Ma = 2.5, T0 = 333 K, R = 3000, alpha = 0.06,
   beta = 0.1 -> omega = 0.0367340 + 0.0005840i (oblique first mode).
   Cases 4 and 5 as tabulated by Sharma & Girimaji (JFM 2023, Table 1);
   Malik used Pr = 0.70, which is what reproduces his values to 7 digits.
6. Ma = 4.5, T0 = 611.11 K, Pr = 0.70, R = 1500, omega = 0.23 (second mode,
   spatial): alpha = 0.2534420 - 0.0027738i (Tumin 2007, JFM 586, quoted by
   Xi, Ren & Fu, arXiv:2006.05970, Table 3, whose own solver gets
   0.253442 - 0.002780i). Malik (1990) lists 0.2534048 - 0.0024921i for
   this case, which neither Tumin, Xi et al. nor this code reproduce.
7. Balakumar & Malik (1992), JFM 239: Ma = 4.5, T0 = 311 K, Pr = 0.72,
   R = 1000, omega = 0.2 -> alpha = 0.220 - 0.003091i (also Tumin 2007).
8. Internal consistency at the DNS conditions (Ma = 0.6, T_inf = 225 K,
   isothermal wall at T_inf, Re_delta* = 1000): the spatial branch fed the
   complex omega of the temporal branch returns alpha, and the Newton
   refinement agrees with the global solve.

In 4-7 the length is the Blasius length sqrt(nu_inf x / U_inf), R = sqrt(Re_x),
the wall is adiabatic for the base flow with T' = 0 for the perturbation, and
viscosity follows Sutherland with S = 110.4 K. The mode is picked by physical
criteria (most unstable, phase speed in a physical range), not by proximity
to the reference.
"""

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")

import numpy as np

from pp.compressibleStability import (
    CompressibleStability,
    Gas,
    compressible_blasius,
    displacement_thickness,
)
from pp.orrSommerfeld import OrrSommerfeld

BLASIUS_FILE = "/home/victor/Desktop/orrSommerfeldSolver/data/blasius.dat"

_fail = 0


def check(name: str, got: complex, want: complex, tol: float) -> None:
    global _fail
    err = abs(got - want)
    ok = err <= tol
    _fail += not ok
    print(
        f"  [{'PASS' if ok else 'FAIL'}] {name}: {complex(got):.7f} "
        f"(reference {complex(want):.7f}, |diff| {err:.1e}, tol {tol:g})"
    )


def most_unstable_temporal(cs: CompressibleStability, alpha: float) -> complex:
    w, _ = cs.temporal(alpha, vectors=False)
    c = w.real / alpha
    # c < 0.9 keeps out the free stream (c = 1) continuous spectrum
    w = w[(c > 0.2) & (c < 0.9)]
    return complex(w[np.argmax(w.imag)])


def most_unstable_spatial(cs: CompressibleStability, omega: float) -> complex:
    a, _ = cs.spatial(omega, vectors=False)
    with np.errstate(divide="ignore", invalid="ignore"):
        c = omega / a.real
    a = a[(a.real > 0) & (c > 0.3) & (c < 1.0) & (np.abs(a.imag) < 0.5 * a.real)]
    return complex(a[np.argmin(a.imag)])


def similarity(mach: float, t0: float, pr: float, twall=None):
    """Adiabatic (twall=None) similarity profile at stagnation temperature t0."""
    gas = Gas(mach=mach, pr=pr, tinf=t0 / (1.0 + 0.2 * mach**2), sutherland=110.4)
    y, u, T = compressible_blasius(gas, twall=twall)
    d99 = float(y[np.argmin(np.abs(u - 0.99))])
    return gas, y, u, T, d99


# --------------------------------------------------------------------------
def operator_vs_sympy() -> None:
    import sympy as sp

    print("1. operator vs symbolic linearisation of compressible NS")
    x, y, z, t, eps = sp.symbols("x y z t epsilon", real=True)
    al, be, om = 0.23 + 0.01j, 0.17, 0.08 - 0.003j
    gam, M, Re, Pr, s = 1.4, 0.8, 700.0, 0.72, 110.0 / 225.0

    Ub = sp.tanh(y)
    Tb = 1 + sp.Rational(3, 4) * sp.exp(-(y**2) / 2)
    qh = [
        sp.exp(-y) * sp.sin(y) * (1 + sp.I * y / 3),
        y * sp.exp(-y),
        sp.exp(-(y**2) / 3) * (sp.I + y),
        sp.exp(-y) * y**2 / (1 + y),
        sp.cos(y) * sp.exp(-y / 2) * (2 - sp.I),
    ]
    E = sp.exp(sp.I * (al * x + be * z - om * t))
    u = Ub + eps * qh[0] * E
    v = eps * qh[1] * E
    w = eps * qh[2] * E
    T = Tb + eps * qh[3] * E
    p = 1 / (gam * M**2) + eps * qh[4] * E
    rho = gam * M**2 * p / T
    mu = T ** sp.Rational(3, 2) * (1 + s) / (T + s)
    lam = -sp.Rational(2, 3) * mu
    X, V = [x, y, z], [u, v, w]
    div = sum(sp.diff(V[i], X[i]) for i in range(3))
    tau = [
        [
            mu * (sp.diff(V[i], X[j]) + sp.diff(V[j], X[i]))
            + (lam * div if i == j else 0)
            for j in range(3)
        ]
        for i in range(3)
    ]

    def D(f):
        return sp.diff(f, t) + sum(V[j] * sp.diff(f, X[j]) for j in range(3))

    eqs = [
        rho * D(V[i])
        + sp.diff(p, X[i])
        - sum(sp.diff(tau[i][j], X[j]) for j in range(3)) / Re
        for i in range(3)
    ]
    Phi = sum(tau[i][j] * sp.diff(V[i], X[j]) for i in range(3) for j in range(3))
    cond = sum(sp.diff(mu * sp.diff(T, X[j]), X[j]) for j in range(3))
    eqs.append(
        rho * D(T)
        + (gam - 1) * rho * T * div
        - gam / (Re * Pr) * cond
        - gam * (gam - 1) * M**2 / Re * Phi
    )
    # continuity, divided by the base density 1 / Tb as in the solver
    eqs.append(
        (sp.diff(rho, t) + sum(sp.diff(rho * V[j], X[j]) for j in range(3))) * Tb
    )
    at0 = {x: 0, z: 0, t: 0}
    lin = [
        sp.lambdify(y, (sp.diff(e, eps).subs(eps, 0) / E).subs(at0), "numpy")
        for e in eqs
    ]
    qf = [sp.lambdify(y, q, "numpy") for q in qh]

    yy = np.linspace(0, 40, 4001)
    gas = Gas(mach=M, gamma=gam, pr=Pr, tinf=225.0, sutherland=110.0)
    cs = CompressibleStability(
        yy, np.tanh(yy), 1 + 0.75 * np.exp(-(yy**2) / 2), Re, gas,
        beta=be, n=161, ymax=40, yi=4, fit_degree=100,
    )
    n = cs.n
    q = np.concatenate([np.asarray(f(cs.yn), dtype=complex) * np.ones(n) for f in qf])
    r = cs.operator(al, om) @ q
    names = ["x momentum", "y momentum", "z momentum", "energy", "continuity"]
    for k in range(5):
        ref = np.asarray(lin[k](cs.yn), dtype=complex) * np.ones(n)
        # boundary rows of the first four equations hold the BCs instead
        sl = slice(1, n - 1) if k < 4 else slice(0, n)
        err = np.max(np.abs(r[k * n : (k + 1) * n][sl] - ref[sl]))
        err /= np.max(np.abs(ref[sl]))
        check(f"{names[k]} relative residual", err, 0.0, 1e-10)


def incompressible_limit() -> None:
    print("2. Ma = 0 against the Orr-Sommerfeld solver (Blasius, Re_delta* = 1000)")
    d = np.loadtxt(BLASIUS_FILE, skiprows=2)
    yb, ub = d[:, 0], d[:, 1]
    os_ = OrrSommerfeld(yb, ub, re=1000.0, n=121)
    cs = CompressibleStability(
        yb, ub, np.ones_like(ub), 1000.0, Gas(mach=0.0, viscosity="constant"), n=121
    )
    w, _ = os_.temporal(0.2)
    w = w[np.abs(w / 0.2) < 1.5]
    check("temporal omega, alpha = 0.2", most_unstable_temporal(cs, 0.2),
          w[np.argmax(w.imag)], 1e-7)
    from pp.localNfactor import select_mode

    a, _ = os_.spatial(0.06, vectors=False)
    check("spatial alpha, omega = 0.06", most_unstable_spatial(cs, 0.06),
          a[select_mode(a, 0.06)], 1e-7)


def jordinson() -> None:
    print("3. Jordinson (1970), spatial Blasius, Re_delta* = 998, omega = 0.1122")
    gas = Gas(mach=0.0)
    y, u, T = compressible_blasius(gas)
    y = y / displacement_thickness(y, u, T)  # delta* units
    cs = CompressibleStability(y, u, T, 998.0, gas, n=121)
    a = most_unstable_spatial(cs, 0.1122)
    # Jordinson quotes four decimals: half a unit in the last place
    check("alpha_r", a.real, 0.3086, 5e-5)
    check("alpha_i", a.imag, -0.0057, 5e-5)


def malik_temporal() -> None:
    print("4. Malik (1990) case 1: Ma = 0.5, R = 2000, alpha = 0.1")
    gas, y, u, T, d99 = similarity(0.5, 278.0, 0.70)
    cs = CompressibleStability(y, u, T, 2000.0, gas, n=141, ymax=25 * d99, fit_degree=70)
    check("omega", most_unstable_temporal(cs, 0.1), 0.0290817 + 0.0022441j, 1e-6)

    print("5. Malik (1990) case 3: Ma = 2.5, R = 3000, alpha = 0.06, beta = 0.1")
    gas, y, u, T, d99 = similarity(2.5, 333.0, 0.70)
    cs = CompressibleStability(
        y, u, T, 3000.0, gas, beta=0.1, n=141, ymax=25 * d99, fit_degree=70
    )
    # Sharma & Girimaji's own solver is 4.8e-6 away from Malik on this case
    check("omega", most_unstable_temporal(cs, 0.06), 0.0367340 + 0.0005840j, 5e-6)


def hypersonic_spatial() -> None:
    print("6. Tumin (2007): Ma = 4.5, T0 = 611.11 K, R = 1500, omega = 0.23")
    gas, y, u, T, d99 = similarity(4.5, 611.11, 0.70)
    cs = CompressibleStability(
        y, u, T, 1500.0, gas, n=161, ymax=3 * d99, yi=0.5 * d99, fit_degree=100
    )
    check("alpha", most_unstable_spatial(cs, 0.23), 0.2534420 - 0.0027738j, 2e-6)

    print("7. Balakumar & Malik (1992): Ma = 4.5, T0 = 311 K, R = 1000, omega = 0.2")
    gas, y, u, T, d99 = similarity(4.5, 311.0, 0.72)
    cs = CompressibleStability(
        y, u, T, 1000.0, gas, n=161, ymax=3 * d99, yi=0.5 * d99, fit_degree=100
    )
    a = most_unstable_spatial(cs, 0.2)
    check("alpha_r (quoted to 3 digits)", a.real, 0.220, 5e-4)
    check("alpha_i", a.imag, -0.003091, 1e-6)


def dns_conditions() -> None:
    print("8. consistency at the DNS conditions (Ma = 0.6, isothermal wall, Re = 1000)")
    gas = Gas(mach=0.6)  # gamma 1.4, Pr 0.72, T_inf = 225 K, S = 110 K
    y, u, T = compressible_blasius(gas, twall=1.0)
    y = y / displacement_thickness(y, u, T)
    cs = CompressibleStability(y, u, T, 1000.0, gas, n=121)
    om = most_unstable_temporal(cs, 0.2)
    a, _ = cs.spatial(om, vectors=False)
    back = complex(a[np.argmin(np.abs(a - 0.2))])
    check("spatial(temporal(alpha = 0.2)) round trip", back, 0.2, 1e-8)
    a0 = most_unstable_spatial(cs, 0.06)
    r = cs.spatial_newton(0.06, a0 * (1.0 + 0.01))
    got = complex(np.nan) if r is None else r[0]
    check("Newton refinement vs global solve", got, a0, 1e-9)


def main() -> None:
    operator_vs_sympy()
    incompressible_limit()
    jordinson()
    malik_temporal()
    hypersonic_spatial()
    dns_conditions()
    print(f"\n{'all checks passed' if _fail == 0 else f'{_fail} check(s) FAILED'}")
    raise SystemExit(1 if _fail else 0)


if __name__ == "__main__":
    main()
