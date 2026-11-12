"""
Spatial local stability (2D TS waves, parallel flow) of the flat-plate boundary layer at the Mach numbers of
the compressible gap DNS, plus the incompressible one, on a (Re_delta*, F) grid.

Thermodynamics as in the Nektar++ sessions (comGapRe1000Ma*/baseflow/dns/session.xml): gamma = 1.4, Pr = 0.72,
Sutherland mu ~ T^1.5 (1 + S) / (T + S) with S = 110 K / Tref, Tref = 225 K (VariableConverter: 110 / m_Tref),
isothermal wall Twall = TInf (T' = 0 at the wall). Base flow: compressible similarity solution with that wall
(pp.compressibleStability.compressible_blasius), lengths in its own displacement thickness delta*.
Incompressible: Ma = 0 with constant viscosity, which is Orr-Sommerfeld exactly (tests/testCompressibleStability.py).

F = omega nu / U^2 * 1e6 = omega / Re_delta* * 1e6 (omega in delta*/U units). For every (Re, F) the TS mode is the
most amplified spatial eigenvalue with phase speed omega / alpha_r in (0.15, 0.75); this excludes the continuous
spectrum (vorticity/entropy c ~ 1, acoustic c = 1 +- 1/Ma).

Output: data/lst/flatPlateMach.npz  (mach, re, F, alpha[mach, re, F])
usage:  python scripts/lstFlatPlateMach.py [--nproc 16]
"""
import argparse
import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

from multiprocessing import Pool

import numpy as np

from pp.compressibleStability import CompressibleStability, Gas, compressible_blasius, displacement_thickness

OUT = os.path.expanduser("~/Desktop/PhD/data/lst/flatPlateMach.npz")
MACH = [0.0, 0.05, 0.2, 0.4, 0.6, 0.8]  # 0.0 = incompressible
RE = np.geomspace(250, 4000, 48)
F = np.linspace(4, 320, 80)
N, YMAX = 81, 60.0
C_RANGE = (0.15, 0.75)


def gas_of(M):
    return Gas(mach=0.0, viscosity="constant") if M == 0.0 else Gas(mach=M)


def profile(M):
    g = gas_of(M)
    y, u, T = compressible_blasius(g, twall=1.0)
    return y / displacement_thickness(y, u, T), u, T


def select(a, om):
    with np.errstate(all="ignore"):
        c = om / a.real
    a = a[(a.real > 0) & (c > C_RANGE[0]) & (c < C_RANGE[1]) & (np.abs(a.imag) < 0.5 * a.real)]
    return complex(a[np.argmin(a.imag)]) if a.size else complex(np.nan, np.nan)


def column(args):
    M, re = args
    y, u, T = profile(M)
    cs = CompressibleStability(y, u, T, re, gas_of(M), n=N, ymax=YMAX)
    out = []
    for f in F:
        om = f * 1e-6 * re
        a, _ = cs.spatial(om, vectors=False)
        out.append(select(a, om))
    return np.array(out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nproc", type=int, default=16)
    a = ap.parse_args()
    jobs = [(M, re) for M in MACH for re in RE]
    with Pool(a.nproc) as p:
        cols = p.map(column, jobs, chunksize=1)
    alpha = np.array(cols).reshape(len(MACH), len(RE), len(F))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    np.savez(OUT, mach=np.array(MACH), re=RE, F=F, alpha=alpha)
    print("saved", OUT, alpha.shape, "nan:", int(np.isnan(alpha).sum()))
