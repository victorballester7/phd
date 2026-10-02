"""
Benchmarks for the Chebyshev Orr-Sommerfeld solver in ``pp.orrSommerfeld``.

Run directly:  python tests/testOrrSommerfeld.py

Each check is a number that does not depend on this code being right:

1. Plane Poiseuille flow at Re = 10000, alpha = 1. Orszag's (1971) eigenvalue
   c = 0.23752649 + 0.00373967i is quoted to eight digits everywhere.
2. The spatial branch fed the *complex* omega that the temporal branch returned
   at alpha = 1 must give alpha = 1 back. This catches any slip in the quartic
   to quadratic reduction or the companion pencil, independently of 1.
3. Blasius: the critical Reynolds number of the parallel flow problem,
   Re_delta* = 519.4, and branch I at Re_delta* = 884.5, which the C++ solver
   in ``orrSommerfeldSolver`` puts at omega = 0.061221, alpha = 0.183792 by a
   temporal sweep plus a second order Gaster transformation.
"""

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
from scipy.optimize import brentq

from pp.orrSommerfeld import OrrSommerfeld

BLASIUS_FILE = "/home/victor/Desktop/orrSommerfeldSolver/data/blasius.dat"
ORSZAG = 0.23752649 + 0.00373967j

_fail = 0


def check(name: str, got: float, want: float, tol: float) -> None:
    global _fail
    ok = abs(got - want) <= tol
    _fail += not ok
    print(
        f"  [{'PASS' if ok else 'FAIL'}] {name}: {got:.8f} "
        f"(reference {want:.8f}, tol {tol:g})"
    )


def poiseuille() -> OrrSommerfeld:
    y = np.linspace(-1.0, 1.0, 2001)
    return OrrSommerfeld(y, 1.0 - y**2, re=10000.0, n=101, ymax=1.0, yi=0.5,
                         fit_degree=40)


def blasius(re: float, n: int = 101) -> OrrSommerfeld:
    d = np.loadtxt(BLASIUS_FILE, skiprows=2)
    return OrrSommerfeld(d[:, 0], d[:, 1], re=re, n=n, ymax=45.0)


def least_damped_temporal(os_: OrrSommerfeld, alpha: float) -> complex:
    w, _ = os_.temporal(alpha)
    w = w[np.abs(w / alpha) < 1.5]
    return complex(w[np.argmax(w.imag)])


def ts_spatial(os_: OrrSommerfeld, omega: float) -> complex | None:
    from pp.localNfactor import select_mode

    al, _ = os_.spatial(omega, vectors=False)
    k = select_mode(al, omega)
    return None if k is None else complex(al[k])


def main() -> None:
    b = blasius(1000.0)
    om = 0.26
    al = ts_spatial(b, om)


if __name__ == "__main__":
    main()
