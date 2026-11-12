"""One impulse run per independent awgn source of the noise runs, to rebuild their exact statistics.

The noise runs force with b_c(x) awgn(1.0) (c = u, v), and Nektar draws awgn independently at
every quadrature point x_p and for each component. In the semi-discrete problem this is

    dq/dt = L q + sum_{p,c} g_{p,c} xi_{p,c}(t),   g_{p,c} = G e_{p,c} b_c(x_p),

with independent white signals xi_{p,c}, so its stationary covariance (the T -> infinity limit of
the time average) is the sum over sources of one impulse integral each:

    C = sigma^2 sum_{p,c} int_0^inf (e^{L t} g_{p,c}) (e^{L t} g_{p,c})^T dt.

Run (p, c) starts from b_c(x_p) at the node x_p only (a window exp(-1e6 r^2), which is 1 at x_p and
0 at every other quadrature point; tested: shifting it by 0.01 gives 1e-44), with no forcing.
Summing the raw second moments <u^2 + v^2> of all runs gives the noise run's infinite-time
statistics, with no statistical error.

Only the element containing (x0, y0) = (-95, 0.5) carries the w = 100 Gaussian; the sources are its
GLL nodes (NUMMODES + 1 per direction), ranked by (w_p b_c(x_p))^2 and kept until 99.98% of that
total (the rest have b < 3e-5 of the largest):
    gap meshes (P5, element [-95.59, -91.21] x [0.25, 0.75]): 5 sources, 88% in v at (-95.21, 0.5)
    flat plate (P6, element [-95.49, -93.24] x [0.278, 0.703]): 10 sources, u and v mixed

Usage: python makeImpulseSourceCases.py   (sessions from <case>_impulse100, see makeImpulseCases.py)
"""

import os
import re

import numpy as np
from numpy.polynomial import legendre as L

SRC = os.path.expanduser("~/Desktop/PhD/src")
GAP = os.path.join(SRC, "incGapRe1000/directLinearSolver/blowingSuction")
FLAT = os.path.join(SRC, "flatPlateRe1000inc/directLinearSolver")

# (directory, impulse100 template, output prefix, element [xa, xb, ya, yb], quadrature points per direction)
CASES = [
    (GAP, "d1.5_w25_impulse100", "d1.5_w25", (-95.5852705, -91.2146882, 0.25, 0.749999985), 7),
    (FLAT, "blowingSuction_impulse100", "blowingSuction", (-95.4918033, -93.2377049, 0.27764403, 0.703152095), 8),
]
X0, Y0, WG = -95.0, 0.5, 100.0
KEEP = 0.9998
B = {  # the Noise function of the noise runs without awgn (A = 0.5/weight)
    "u": "-2*{w}*(y-({y0}))*{a}*exp(-{w}*((x-({x0}))^2 + (y-({y0}))^2))",
    "v": "2*{w}*(x-({x0}))*{a}*exp(-{w}*((x-({x0}))^2 + (y-({y0}))^2))",
}


def gll(n):
    c = np.zeros(n)
    c[-1] = 1
    x = np.concatenate(([-1.0], np.sort(L.legroots(L.legder(c))), [1.0]))
    return x, 2 / (n * (n - 1) * L.legval(x, c) ** 2)


def sources(elem, nq):
    xa, xb, ya, yb = elem
    xi, wi = gll(nq)
    X, Y = (xa + xb) / 2 + xi * (xb - xa) / 2, (ya + yb) / 2 + xi * (yb - ya) / 2
    W = np.outer(wi * (xb - xa) / 2, wi * (yb - ya) / 2)
    A = 0.5 / WG
    out = []
    for i in range(nq):
        for j in range(nq):
            e = np.exp(-WG * ((X[i] - X0) ** 2 + (Y[j] - Y0) ** 2))
            for c, b in (("u", -2 * WG * (Y[j] - Y0) * A * e), ("v", 2 * WG * (X[i] - X0) * A * e)):
                out.append(((W[i, j] * b) ** 2, c, X[i], Y[j], b))
    out.sort(reverse=True)
    tot = sum(o[0] for o in out)
    kept, cum = [], 0.0
    for o in out:
        if cum >= KEEP:
            break
        cum += o[0] / tot
        kept.append((o[1], o[2], o[3], o[4], o[0] / tot))
    return kept


def make(base, template, prefix, elem, nq):
    s0 = open(os.path.join(base, template, "session.xml")).read()
    for k, (c, xp, yp, b, share) in enumerate(sources(elem, nq)):
        name = f"{prefix}_src{k}{c}"
        win = f"exp(-1e6*((x-({xp:.9f}))^2 + (y-({yp:.9f}))^2))"
        shape = B[c].format(w=f"{WG:g}", a=f"{0.5 / WG:g}", x0=f"{X0:g}", y0=f"{Y0:g}")
        ic = {v: "0" for v in "uvp"}
        ic[c] = f"{shape}*{win}"
        new = (
            f"<!-- Impulse of one awgn source of the noise runs: component {c} at the quadrature node\n"
            f"         ({xp:.7f}, {yp:.7f}), value b_{c} = {b:.4e} ({100 * share:.2f}% of the forcing variance);\n"
            "         the window is 1 at that node and 0 at every other quadrature point (makeImpulseSourceCases.py) -->\n"
            '    <FUNCTION NAME="InitialConditions">\n'
            + "".join(f'      <E VAR="{v}" VALUE="{ic[v]}" />\n' for v in "uvp")
            + "    </FUNCTION>"
        )
        s, n = re.subn(r'(?<!<!-- )<FUNCTION NAME="InitialConditions">\s*<E VAR="u" VALUE="-2\*.*?</FUNCTION>',
                       new, s0, count=1, flags=re.S)
        assert n == 1, name
        os.makedirs(os.path.join(base, name), exist_ok=True)
        open(os.path.join(base, name, "session.xml"), "w").write(s)
        print(f"{name}: {c} at ({xp:.7f}, {yp:.7f}), b = {b: .3e}, share {100 * share:.2f}%")


if __name__ == "__main__":
    for case in CASES:
        make(*case)
