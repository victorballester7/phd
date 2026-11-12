"""Create restarts of the quasi-3D d=1.5 runs with a controlled spanwise-varying seed.

The unforced runs only carry the 3D disturbance left over from the initial noise, which
decays to round-off, so the breakdown position is set by that numerical floor. These
cases add a small body force in w, localised just downstream of the gap and containing
no spanwise-uniform part, so that the 2D limit cycle is left untouched while the 3D
seed level is fixed by the amplitude A:

    f_w = A exp(-((x - xF)^2/16 + (y - 1)^2/0.5)) Z(z) T(t)

Z is the sum of the first four spanwise harmonics of Lz and T a multisine of 30 lines,
omega = 0.01..0.30 (seed 1), both normalised to unit rms. If breakdown is seed limited,
x_b moves upstream as A grows; if not, it stays at x - w ~ 450-490.

Usage: python makeSeedForcingCases.py   (edit CASES below)
"""

import os
import re
import shutil
import numpy as np

SRC = os.path.expanduser("~/Desktop/PhD/src/quasi3DGapRe1000inc/baseflow/dns")

# (parent case, width, restart checkpoint number, amplitudes)
# 2026-10-26: w70 from chk 200 and w66 from chk 275, amplitudes 1e-5 and 1e-3 (1e-9 added 10-29).
# 2026-11-08: region III, to test the local secondary stability (docs/notes/20261108_secondaryStability.md):
#             w63 from chk 39 and w64 from chk 59, amplitudes 1e-5 and 1e-9.
# 2026-11-11: w60 from chk 59 (t = 5507, the last checkpoint of the unforced run), 1e-5 and 1e-9.
CASES = [
    ("d1.5_w60_Lz6.0_k32", 60, 59, ["1e-5", "1e-9"]),
]
X_OFFSET = 20.0  # forcing centre, distance downstream of the gap trailing edge
N_LINES = 30
D_OMEGA = 0.01
N_HARMONICS = 4


def forcing_expression():
    rng = np.random.default_rng(1)
    psi = rng.uniform(0.0, 2.0 * np.pi, N_LINES)
    phi = rng.uniform(0.0, 2.0 * np.pi, N_HARMONICS)
    T = " + ".join(f"sin({D_OMEGA * (m + 1):.2f}*t + {psi[m]:.4f})" for m in range(N_LINES))
    Z = " + ".join(f"cos({k + 1}*2*PI/LZ*z + {phi[k]:.4f})" for k in range(N_HARMONICS))
    # unit rms: a sum of n unit sines has variance n/2
    return (
        f"seedAmp*exp(-((x-xForce)^2/16 + (y-1)^2/0.5))"
        f" * ({Z})/sqrt({N_HARMONICS / 2})"
        f" * ({T})/sqrt({N_LINES / 2})"
    )


def make_case(parent, width, nchk, amp):
    name = f"{parent}_seed{amp}"
    dst = os.path.join(SRC, name)
    os.makedirs(dst, exist_ok=True)
    for f in ("mesh.geo", "mesh.msh", "mesh.xml", "pbspro.job", "slurm.job"):
        shutil.copy(os.path.join(SRC, parent, f), dst)

    with open(os.path.join(SRC, parent, "session.xml")) as fh:
        s = fh.read()

    s = re.sub(r'<F VAR="u,v,w,p" FILE="mesh_\d+\.chk" />',
               f'<F VAR="u,v,w,p" FILE="mesh_{nchk}.chk" />', s, count=1)
    params = (
        "      <!-- spanwise-varying seed, see scripts/makeSeedForcingCases.py -->\n"
        f"      <P> seedAmp        = {amp}                </P>\n"
        f"      <P> xForce         = {width + X_OFFSET}   </P>\n"
    )
    s = s.replace("      <P> HomModesZ      = 32                 </P>\n",
                  "      <P> HomModesZ      = 32                 </P>\n" + params, 1)
    func = (
        '<FUNCTION NAME="BodyForce">\n'
        f'  <E VAR="w" VALUE="{forcing_expression()}" />\n'
        "</FUNCTION>\n"
    )
    s = s.replace("  </CONDITIONS>", func + "  </CONDITIONS>\n"
                  "  <FORCING>\n"
                  '    <FORCE TYPE="Body">\n'
                  "      <BODYFORCE> BodyForce </BODYFORCE>\n"
                  "    </FORCE>\n"
                  "  </FORCING>", 1)
    assert "seedAmp" in s and "<FORCING>" in s and f"mesh_{nchk}.chk" in s
    with open(os.path.join(dst, "session.xml"), "w") as fh:
        fh.write(s)

    tag = f"w{width}s{amp.replace('e-', 'm')}"
    for job in ("pbspro.job", "slurm.job"):
        p = os.path.join(dst, job)
        with open(p) as fh:
            j = fh.read()
        j = re.sub(r"(#PBS -N |--job-name=)\S+", rf"\g<1>{tag}", j)
        with open(p, "w") as fh:
            fh.write(j)
    print(f"created {dst}")


if __name__ == "__main__":
    for parent, width, nchk, amps in CASES:
        for amp in amps:
            make_case(parent, width, nchk, amp)
