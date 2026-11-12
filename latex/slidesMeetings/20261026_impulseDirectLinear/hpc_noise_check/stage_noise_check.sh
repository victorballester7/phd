#!/bin/bash
# Rank-1 white-noise run on the gap (Ma 0.05, d1.5_w25, mesh R1), to compare with the impulse run meshtest_R1.
# Run on the HPC login node from ~/nektar-jfdirect: ./stage_noise_check.sh <rebuild job id>
set -e
R=$HOME/Desktop/PhD/runs/comGapRe1000Ma0.05/directLinearSolver/blowingSuction/d1.5_w25
I=$R/meshtest_R1; D=$R/meshtest_R1_noise
mkdir -p $D && cd $D
ln -sf $I/mesh.xml mesh.xml
ln -sfn $I/baseflow.fld baseflow.fld
cp -r $I/baseflow_evolved.fld . 2>/dev/null || true   # same T, dt, base flow: the cached base evolution is reused
python3 - $I/session.xml session.xml <<'PY'
import re, sys
s = open(sys.argv[1]).read()
def sub1(a, b, flags=0):
    global s
    s2, n = re.subn(a, b, s, count=1, flags=flags); assert n == 1, a; s = s2
# t in [0, 23660]: transient up to 3640 (one transit, as the impulse run), statistics over [3640, 23660]
sub1(r'<P> NumSteps .*?</P>', '<P> NumSteps        = 118300     </P> <!-- final time 23660 -->')
sub1(r'<P> JacobiFreeStatsStartTime .*?</P>', '<P> JacobiFreeStatsStartTime  = 3640 </P> <!-- averages over [3640, 23660] -->')
# the vortex b becomes the forcing shape; zero initial perturbation
blk = re.search(r'^[ \t]*<FUNCTION NAME="InitialConditions">.*?</FUNCTION>', s, flags=re.S | re.M).group(0)  # not the commented one
assert 'exp(-weight' in blk
dist = blk.replace('NAME="InitialConditions"', 'NAME="Disturbance"')
zero = '''<FUNCTION NAME="InitialConditions">
      <E VAR="rho"  VALUE="0" />
      <E VAR="rhou" VALUE="0" />
      <E VAR="rhov" VALUE="0" />
      <E VAR="E"    VALUE="0" />
    </FUNCTION>'''
s = s.replace(blk, zero + '\n\n    <!-- forcing shape b: the vortex of the impulse run -->\n    ' + dist)
# rank-1 white noise: one random number per time level for the whole field (FIELDFORCETIMEFCN), amplitude
# a = 1/sqrt(0.586 dt) so that sigma^2 = a^2 dt (b1^2 + b2^2) = 1 (DIRK2 samples it at its two stage times)
sub1(r'(\s*)</FORCING>', r'''\1  <!-- rank-1 white noise: f = a xi(t) b, xi one N(0,1) value per time level; sigma^2 = 1 -->
    <FORCE TYPE="Field" PERTURBATION="True">
      <FIELDFORCE> Disturbance </FIELDFORCE>
      <FIELDFORCETIMEFCN> awgn(1.0)/sqrt(0.5857864376269049*0.2) </FIELDFORCETIMEFCN>
    </FORCE>\1</FORCING>''')
open(sys.argv[2], 'w').write(s)
PY
sed -e 's/#PBS -N meshR1/#PBS -N noiseR1/' -e 's/#PBS -l walltime=12:00:00/#PBS -l walltime=24:00:00/' $I/pbspro.job > pbspro.job
grep -c "FIELDFORCETIMEFCN" session.xml
qsub ${1:+-W depend=afterok:$1} pbspro.job
