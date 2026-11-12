"""List the compressible gap DNS runs that reached a steady base flow, with their last checkpoint.

A run counts as steady when, at every history point, the velocity (rhou/rho, rhov/rho) changes by less
than TOL over the last WINDOW time units of HistoryPoints.his (.his.gz + .his, see hisparts). Only plain d<D>_w<W> case folders are
considered (variants such as _old, _new, _smallDt are skipped). Pure Python (the HPC login node has no
numpy). Output, one line per steady case: "<Ma> <case> <checkpoint>", the format read by
runs/pp_submit_impulse.sh and runs/pp_prep_case.sh.

usage: python3 steadyBaseflows.py [runs dir] > pp_steady_cases.txt   (unsteady/unclear cases go to stderr)
"""
import gzip
import os
import re
import shlex
import subprocess
import sys

RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Desktop/PhD/runs")
MACH = ["0.05", "0.2", "0.4", "0.6", "0.8"]
TOL = 1e-6
WINDOW = 150.0


def hisparts(his):
    """HistoryPoints.his of a case: finished runs in his.gz (pbspro.sh since 2026-11-10), old scheme his.old, and the
    run in progress in his itself; oldest first."""
    return [p for p in (his + ".gz", his + ".old", his) if os.path.isfile(p)]


def histail(his, nlines):
    """last nlines of the whole history (his.gz + his.old + his) as a list of lines"""
    cmds = [("zcat " if p.endswith(".gz") else "cat ") + shlex.quote(p) for p in hisparts(his)]
    sh = "{ %s; } | tail -n %d" % ("; ".join(cmds), nlines)
    return subprocess.check_output(sh, shell=True).decode().split("\n")


def hisheader(his):
    """leading # lines of the history"""
    p = hisparts(his)[0]
    out = []
    with (gzip.open(p, "rt") if p.endswith(".gz") else open(p)) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            out.append(line)
    return out


def npoints(his):
    return sum(1 for line in hisheader(his) if re.match(r"#\s+\d+\s", line))


def change(his, npts):
    """max |du|, |dv| over the last WINDOW time units, and the time span actually covered"""
    nlines = int(npts * (WINDOW / 0.5 + 50))  # samples every 0.5 or more; take a margin
    tail = histail(his, nlines)
    rows = [list(map(float, l.split())) for l in tail if l and not l.startswith("#")]
    rows = rows[len(rows) % npts:]  # whole samples only
    samples = [rows[i:i + npts] for i in range(0, len(rows), npts)]
    if len(samples) < 2:
        return None, 0.0
    t_end = samples[-1][0][0]
    samples = [s for s in samples if s[0][0] >= t_end - WINDOW]
    last = samples[-1]
    dmax = 0.0
    for s in samples:
        for p in range(npts):
            for k in (2, 3):
                dmax = max(dmax, abs(s[p][k] / s[p][1] - last[p][k] / last[p][1]))
    return dmax, t_end - samples[0][0][0]


def last_chk(case):
    ks = []
    for d in os.listdir(case):
        m = re.fullmatch(r"mesh_(\d+)\.chk", d)
        if m and os.path.exists(os.path.join(case, d, "Info.xml")):
            ks.append(int(m.group(1)))
    return f"mesh_{max(ks)}.chk" if ks else None


for M in MACH:
    base = os.path.join(RUNS, f"comGapRe1000Ma{M}", "baseflow", "dns")
    for c in sorted(os.listdir(base)):
        case = os.path.join(base, c)
        his = os.path.join(case, "HistoryPoints.his")
        if not re.fullmatch(r"d[0-9.]+_w[0-9.]+", c) or not hisparts(his):
            continue
        try:
            dmax, span = change(his, npoints(his))
        except Exception as e:
            print(f"{M} {c}: unreadable history ({e})", file=sys.stderr)
            continue
        chk = last_chk(case)
        if dmax is not None and span >= WINDOW - 1 and dmax < TOL and chk:
            print(f"{M} {c} {chk}")
        else:
            print(f"{M} {c}: not steady (max change {dmax} over {span:.0f} time units, last chk {chk})", file=sys.stderr)
