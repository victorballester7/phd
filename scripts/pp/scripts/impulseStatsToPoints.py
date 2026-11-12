"""Interpolate the statistics of the impulse-response runs onto the noise runs' points.

An impulse run (``scripts/makeImpulseCases.py``) writes cumulative ReynoldsStresses averages
from t = 0: ``mesh_<k>_stress.fld`` every 192 time units and ``mesh_stress.fld`` at the end.
Each is proportional to int_0^t q q^T dt, so N(x) from the dump at time t is N(x) of the
integral up to t, and the change between dumps measures how converged the integral is.

For every case this fetches the dumps from typhoon (rsync, not the sshfs mount) into a
temporary directory, interpolates them locally with FieldConvert onto the parent noise
run's ``data/pointsavg_n600.pts`` and writes

    <impulse dir>/data/pointsavg_n600.dat   final dump (same format as the noise runs)
    <impulse dir>/data/convergence.npz      t, x and the L2 amplitude A(x) of every dump
                                            (the interpolated dumps themselves are not kept)

Usage: python scripts/impulseStatsToPoints.py [case ...]   (default: all impulse dirs found)
"""

import contextlib
import io
import os
import shutil
import re
import subprocess
import sys
import tempfile

import numpy as np

from pp.DeltaN_computation import computeAmplitude

SRC = os.path.expanduser("~/Desktop/PhD/src")
GAP = os.path.join(SRC, "incGapRe1000/directLinearSolver/blowingSuction")
FLAT = os.path.join(SRC, "flatPlateRe1000inc/directLinearSolver")
REMOTE = "typhoon:Desktop/PhD/runs"  # the impulse batch runs on typhoon (2026-10-29)
N = 600


def impulse_dirs():
    out = []
    for base in (GAP, FLAT):
        for d in sorted(os.listdir(base)):
            if re.search(r"(_impulse(100(src)?)?|_src\d+[uv])$", d):
                out.append(os.path.join(base, d))
    return out


def parent_of(path):
    base, name = os.path.split(path)
    return os.path.join(base, "blowingSuction" if base == FLAT else re.sub(r"(_impulse(100(src)?)?|_src\d+[uv])$", "", name))


def final_time(fld):
    info = open(os.path.join(fld, "Info.xml")).read()
    return float(re.search(r"<FinalTime>([^<]*)", info).group(1))


def interp(mesh, fld, pts, out):
    subprocess.run(
        ["FieldConvert", "-f", "-m", f"interppoints:fromxml={mesh}:fromfld={fld}:topts={pts}", out],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )


def process(path):
    name = os.path.basename(path)
    parent = parent_of(path)
    # the _impulse100src runs have their own source-refined mesh (scripts/refineSourceMesh.py); the others use the parent's
    mesh = os.path.join(path, "mesh.xml") if os.path.isfile(os.path.join(path, "mesh.xml")) else os.path.join(parent, "mesh.xml")
    pts = os.path.join(parent, "data", f"pointsavg_n{N}.pts")
    assert os.path.isfile(mesh) and os.path.isfile(pts), (mesh, pts)
    rel = os.path.relpath(path, SRC)
    os.makedirs(os.path.join(path, "data"), exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run(
            ["rsync", "-a", "--include=mesh_stress.fld/***", "--include=mesh_*_stress.fld/***",
             "--include=output.txt", "--exclude=*", f"{REMOTE}/{rel}/", tmp + "/"],
            capture_output=True, text=True,
        )
        if r.returncode != 0 or not os.path.isdir(os.path.join(tmp, "mesh_stress.fld")):
            print(f"{name}: no final mesh_stress.fld yet")
            return
        if os.path.isfile(os.path.join(tmp, "output.txt")):
            shutil.copy(os.path.join(tmp, "output.txt"), os.path.join(path, "data", "output.txt"))
        dumps = [d for d in os.listdir(tmp) if re.fullmatch(r"mesh_\d+_stress\.fld", d)]
        final = os.path.join(path, "data", f"pointsavg_n{N}.dat")
        interp(mesh, os.path.join(tmp, "mesh_stress.fld"), pts, final)
        ts, As = [], []
        for d in dumps:
            dat = os.path.join(tmp, f"{name}.dat")  # the case name keeps extract_depth_width happy
            interp(mesh, os.path.join(tmp, d), pts, dat)
            with contextlib.redirect_stdout(io.StringIO()):
                x, A = computeAmplitude(dat, False, field="rms")
            ts.append(final_time(os.path.join(tmp, d)))
            As.append(A)
        T = final_time(os.path.join(tmp, "mesh_stress.fld"))
        with contextlib.redirect_stdout(io.StringIO()):
            x, A = computeAmplitude(final, False, field="rms")
        ts.append(T)
        As.append(A)
        order = np.argsort(ts)
        np.savez(os.path.join(path, "data", "convergence.npz"), t=np.array(ts)[order], x=x, A=np.array(As)[order])
        print(f"{name}: final t = {T:.0f}, {len(dumps)} intermediate dumps")


if __name__ == "__main__":
    wanted = sys.argv[1:]
    for p in impulse_dirs():
        if not wanted or os.path.basename(p) in wanted:
            process(p)
