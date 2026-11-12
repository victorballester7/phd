"""Create impulse-response counterparts of the incompressible blowingSuction (awgn) runs.

The noise runs force the linearised equations with f(x, t) = b(x) awgn(t), b the Gaussian
vortex of the "Noise" function, and time-average the Reynolds stresses over [5040, 24000].
For a forcing b xi(t) with xi white in time, the stationary covariance is

    <q q^T> = sigma^2 int_0^inf q_imp(t) q_imp(t)^T dt,    q_imp(0) = b, no forcing,

so one impulse run gives the converged statistics: the ReynoldsStresses filter averages from
t = 0, i.e. it accumulates int_0^T q q^T dt / T, and N(x) = log(A/A0) does not depend on 1/T.
It must run until the response has left the domain (or decayed, inside the gap).

Two shapes b per case:
    <case>_impulse      weight 1    resolved vortex (radius 0.7), as d1.5_w25_impulse and the
                                    compressible batch
    <case>_impulse100   weight 100  the noise runs' envelope. On the P5 meshes the element at
                                    x0 = -95 spans [-95.59, -91.21] x [0.25, 0.75], so the
                                    Gaussian (radius 0.07) is only seen by the quadrature column
                                    x = -95.23: a wall-normal point force at y ~ 0.5, which is
                                    also what the per-point awgn forcing reduces to

Usage: python makeImpulseCases.py   (edit CASES below; then rsync to the HPC, link mesh.xml and
baseflow.fld to the parent case and submit)
"""

import os
import re

SRC = os.path.expanduser("~/Desktop/PhD/src")
GAP = os.path.join(SRC, "incGapRe1000/directLinearSolver/blowingSuction")
FLAT = os.path.join(SRC, "flatPlateRe1000inc/directLinearSolver")

# (directory holding the noise run, parent case, output prefix)
CASES = [(GAP, c, c) for c in ("d0.5_w10", "d0.5_w90", "d1.25_w90", "d1.5_w25", "d2_w28", "d2.5_w15", "d4_w15")]
CASES += [(FLAT, "blowingSuction_envL", "blowingSuction")]  # envL: flat-plate session with the IMEX block

WEIGHTS = {"impulse": 1, "impulse100": 100}
SKIP = set()  # d1.5_w25_impulse was first run on the HPC with nektar-jfdirect (t = 3640); rerun with the batch

T_END = 4800  # final time; the response has left x = 975 by t ~ 3300 for d1.5_w25
T_STATS = 192  # cumulative statistics dumps, to check the convergence of the time integral
T_CHK = 960
T_INFO = 12
T_HIS = 1.2


def steps(t, dt):
    n = round(t / dt)
    assert abs(n * dt - t) < 1e-9, (t, dt)
    return n


def make_case(base, parent, prefix, suffix, weight):
    name = f"{prefix}_{suffix}"
    dst = os.path.join(base, name)
    with open(os.path.join(base, parent, "session.xml")) as fh:
        s = fh.read()

    def sub1(pat, new, flags=0):
        nonlocal s
        s, n = re.subn(pat, new, s, count=1, flags=flags)
        assert n == 1, (name, pat)

    dt = float(re.search(r"<P> TimeStep\s*=\s*([0-9.eE+-]+)\s*</P>", s).group(1))
    sub1(r"<P> NumSteps .*?</P>", f"<P> NumSteps        = {steps(T_END, dt)}     </P> <!-- final time {T_END} -->")
    sub1(r"<P> IO_CheckSteps .*?</P>", f"<P> IO_CheckSteps   = {steps(T_CHK, dt)}     </P> <!-- every {T_CHK} -->")
    sub1(r"<P> IO_InfoSteps .*?</P>", f"<P> IO_InfoSteps    = {steps(T_INFO, dt)}        </P>")
    sub1(r"<P> HistorySteps .*?</P>", f"<P> HistorySteps    = {steps(T_HIS, dt)} </P>")
    sub1(r"<P> AverageSteps .*?</P>",
         f"<P> AverageSteps    = {steps(T_STATS, dt)} </P> <!-- cumulative statistics every {T_STATS} -->")

    note = "resolved vortex (radius 0.7)" if weight == 1 else "envelope of the noise runs"
    if re.search(r"<P> weight\s*=", s):
        sub1(r"<P> weight\s*=\s*[0-9.]+\s*</P>",
             f"<P> weight          = {weight}                  </P> <!-- {note} -->")
    else:
        for w in ("weightX", "weightY"):
            sub1(rf"<P> {w}\s*=\s*[0-9.]+\s*</P>", f"<P> {w}         = {weight}                  </P>")
        s = s.replace("<!-- streamwise weight", f"<!-- {note}; streamwise weight", 1)

    # initial condition = forcing shape b: the Noise expressions without the awgn factor
    noise = re.search(r'<FUNCTION NAME="Noise">(.*?)</FUNCTION>', s, flags=re.S).group(1)
    shape = {v: e.replace("*awgn(1.0)", "") for v, e in re.findall(r'<E VAR="(\w)" VALUE="([^"]*)"', noise)}
    assert set(shape) == {"u", "v", "p"} and "awgn" not in "".join(shape.values()), name
    ic = (
        "<!-- Impulse response: the initial perturbation is the forcing shape b (the Noise function without\n"
        "         awgn), with no forcing. For forcing b xi(t), xi white noise of intensity sigma^2,\n"
        "         <q q^T> = sigma^2 int_0^inf q_imp q_imp^T dt = sigma^2 T x (ReynoldsStresses average from t = 0),\n"
        "         once the response has left the domain -->\n"
        '    <FUNCTION NAME="InitialConditions">\n'
        + "".join(f'      <E VAR="{v}" VALUE="{shape[v]}" />\n' for v in ("u", "v", "p"))
        + "    </FUNCTION>"
    )
    # the active (uncommented) InitialConditions block
    sub1(r'(?<!<!-- )<FUNCTION NAME="InitialConditions">\s*(<E VAR="\w" VALUE="0" />\s*){3}</FUNCTION>', ic)
    sub1(r"\s*<FORCING>.*?</FORCING>", "", flags=re.S)
    assert "<FORCING>" not in s and "FIELDFORCE" not in s, name

    os.makedirs(dst, exist_ok=True)
    with open(os.path.join(dst, "session.xml"), "w") as fh:
        fh.write(s)

    tag = re.sub(r"[^0-9a-z]", "", prefix.replace("blowingSuction", "flat"))[:10]
    walltime = "16:00:00" if base == FLAT else "06:00:00"
    with open(os.path.join(dst, "pbspro.job"), "w") as fh:
        fh.write(
            "#!/bin/bash\n"
            "#PBS -l select=1:ncpus=32:mpiprocs=32:mem=32gb\n"
            f"#PBS -l walltime={walltime}\n"
            f"#PBS -N {'i' if weight == 1 else 'J'}{tag}\n"
            "#PBS -e log.job\n"
            "#PBS -o output.job\n"
            "\n"
            "source $HOME/Desktop/PhD/scripts/pbspro.sh\n"
        )
    print(f"{name}: dt = {dt}, NumSteps = {steps(T_END, dt)}, weight = {weight}")
    return name


if __name__ == "__main__":
    for base, parent, prefix in CASES:
        for suffix, weight in WEIGHTS.items():
            if f"{prefix}_{suffix}" not in SKIP:
                make_case(base, parent, prefix, suffix, weight)
