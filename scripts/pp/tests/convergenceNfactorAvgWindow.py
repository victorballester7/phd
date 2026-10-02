#!/usr/bin/env python3
"""Convergence study of the time-averaging window of the stochastic-forcing runs.

The amplification (N factor) curves are built from Reynolds-stress fields that
Nektar accumulates as a *cumulative* time average from t = 0. The production
pipeline (``scripts/combineAvg.py`` -> ``scripts/createPointsOfSectionDomain*``
-> ``pp.DeltaN_computation``) keeps a single window: it subtracts the cumulative
average at the first dump past ``mintime`` (transient removal) from the last
dump, and works with that one field.

This script does the same thing ``NSLICES`` times, keeping the *start* of the
window fixed at ``mintime`` and sliding the *end* over equidistant final times,
so every slice is an average over

    (t_beg, t_beg + k/K * (t_last - t_beg)],   k = 1 .. K

with the k = K slice being exactly the production window. The N factor curves of
all slices are then overlaid with a linear colour map in the number of time
steps used, which is the point of the exercise: if the short-window curves
converge onto the long-window one, the average is statistically converged.

Everything it writes is tagged with the ``convAvg_`` prefix (remote .fld files)
and lives under ``data/convergence/`` locally, so it can never be confused with
the production ``mesh_avg.fld`` / ``data/pointsavg_n*.dat``.

Usage:
    python tests/convergenceAvgWindow.py                       # default case
    python tests/convergenceAvgWindow.py --dir <case dir> -k 10
    python tests/convergenceAvgWindow.py --force               # redo remote work
    python tests/convergenceAvgWindow.py --no-remote           # replot only
"""

import argparse
import json
import os
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from createPointsOfSectionDomainOneChkFile import gen_points_file  # noqa: E402
from pp.colors import colors  # noqa: E402
from pp.DeltaN_computation import computeNx  # noqa: E402
from pp.fileManagement import extract_depth_width, extractValueXML  # noqa: E402
from pp.ssh.utilities import get_remote_dir, run_cmd  # noqa: E402

DEFAULT_DIR = (
    "/home/victor/Desktop/PhD/src/incGapRe1000/directLinearSolver/"
    "blowingSuction/d1.75_w40_old"
)

# Prefix marking every artefact of this study, so it is obvious at a glance
# (locally and on the HPC) that these files are throw-away convergence probes.
CONV_PREFIX = "convAvg"
CONV_SUBDIR = "data/convergence"

STRESS_SUFFIX = "_stress.fld"
MINTIME = 5000.0  # transient cut, same value as scripts/combineAvg.py
NSLICES = 9
NPOINTS_Y = 600  # wall-normal points, same as the production pts file
DO_LOO = False  # False -> L2 amplitude, as in the production N factor curves


# ---------------------------------------------------------------------------
# remote / mounted helpers
# ---------------------------------------------------------------------------
def get_mounted_dir(dir_local: str, host: str = "hpc", user: str = "vb824") -> str:
    """Mounted-filesystem view of the remote case directory (same mapping as
    ``pp.ssh.fieldconvert.combineAvg``), used for cheap file inspection so that
    only the actual FieldConvert calls need SSH."""
    dir_remote = get_remote_dir(dir_local)
    dir_mounted = dir_remote.replace(
        f"/rds/general/user/{user}/home/Desktop", f"/home/victor/hosts/{host}"
    )
    if not os.path.isdir(dir_mounted):
        raise FileNotFoundError(
            colors.FAIL + f"Mounted directory not found: {dir_mounted}" + colors.ENDC
        )
    return dir_mounted


def read_stress_info(folder: str) -> Tuple[float, int]:
    """(FinalTime, |NumberOfFieldDumps|) of one ``*_stress.fld`` dump.

    The absolute value matters: ``combineAvg`` flips the sign of the start
    folder's dump count in place, so a folder that has already served as a
    window start on a previous run carries a negative count.
    """
    root = ET.parse(os.path.join(folder, "Info.xml")).getroot()
    final_time = float(root.findtext(".//FinalTime", default="0"))
    dumps = abs(int(root.findtext(".//NumberOfFieldDumps", default="0")))
    return final_time, dumps


def scan_stress_folders(dir_mounted: str) -> List[Tuple[str, float, int]]:
    """All stress dumps as (name, FinalTime, dumps), sorted by FinalTime."""
    out = []
    for d in sorted(os.listdir(dir_mounted)):
        if not d.endswith(STRESS_SUFFIX):
            continue
        path = os.path.join(dir_mounted, d)
        if not os.path.exists(os.path.join(path, "P0000000.fld")):
            continue
        try:
            final_time, dumps = read_stress_info(path)
        except Exception as e:
            print(colors.WARNING + f"Skipping {d}: {e}" + colors.ENDC)
            continue
        out.append((d, final_time, dumps))
    out.sort(key=lambda t: t[1])
    return out


def patch_start_folder(dir_mounted: str, folder: str) -> None:
    """Negate NumberOfFieldDumps of the window-start dump, which is how
    FieldConvert's combineAvg is told to *subtract* it. Idempotent: an already
    negative value is left alone."""
    xmlfile = os.path.join(dir_mounted, folder, "Info.xml")
    tree = ET.parse(xmlfile)
    node = tree.getroot().find(".//NumberOfFieldDumps")
    if node is None or node.text is None:
        print(colors.WARNING + f"No NumberOfFieldDumps in {xmlfile}" + colors.ENDC)
        return
    if node.text.strip().startswith("-"):
        print(
            colors.OKBLUE
            + f"Start folder {folder} already patched ({node.text})"
            + colors.ENDC
        )
        return
    node.text = str(-int(node.text))
    tree.write(xmlfile)
    print(colors.WARNING + f"Patched NumberOfFieldDumps -> {node.text}" + colors.ENDC)


def read_run_parameters(dir_local: str) -> Tuple[float, int]:
    """(TimeStep, SampleFrequency) of the run, from session.xml."""
    session = os.path.join(dir_local, "session.xml")
    dt = float(extractValueXML(session, "TimeStep"))
    with open(session) as f:
        m = re.search(
            r'<PARAM\s+NAME\s*=\s*"SampleFrequency"\s*>\s*([0-9]+)', f.read()
        )
    sample_frequency = int(m.group(1)) if m else 1
    if m is None:
        print(
            colors.WARNING
            + "SampleFrequency not found in session.xml, assuming 1"
            + colors.ENDC
        )
    return dt, sample_frequency


# ---------------------------------------------------------------------------
# slice selection
# ---------------------------------------------------------------------------
class Slice:
    """One averaging window (t_beg, t_end] and everything derived from it."""

    def __init__(self, k: int, folder: str, t_beg: float, t_end: float,
                 dumps_beg: int, dumps_end: int, dt: float, sample_frequency: int):
        self.k = k
        self.folder = folder
        self.t_beg = t_beg
        self.t_end = t_end
        self.nsamples = dumps_end - dumps_beg
        self.nsteps = self.nsamples * sample_frequency
        self.T_avg = t_end - t_beg
        self.fld = f"{CONV_PREFIX}_k{k:02d}.fld"
        self.datbase = f"{CONV_SUBDIR}/pointsavg_n{NPOINTS_Y}_{CONV_PREFIX}_k{k:02d}"

    def __repr__(self) -> str:
        return (
            f"k={self.k:2d}  {self.folder:<24s}  T=({self.t_beg:8.1f}, {self.t_end:8.1f}]"
            f"  T_avg={self.T_avg:8.1f}  samples={self.nsamples:8d}  steps={self.nsteps:9d}"
        )


def build_slices(
    stress: List[Tuple[str, float, int]],
    mintime: float,
    nslices: int,
    dt: float,
    sample_frequency: int,
) -> Tuple[str, List[Slice]]:
    """Fixed window start at the first dump past ``mintime``; ``nslices``
    equidistant window ends spanning up to the last available dump."""
    beg_idx = next((i for i, (_, t, _) in enumerate(stress) if t >= mintime), None)
    if beg_idx is None:
        raise RuntimeError(
            colors.FAIL + f"No stress folder with FinalTime >= {mintime}" + colors.ENDC
        )
    beg_folder, t_beg, dumps_beg = stress[beg_idx]
    print(
        colors.OKGREEN
        + f"Window start: {beg_folder} (FinalTime = {t_beg:.1f}, dumps = {dumps_beg})"
        + colors.ENDC
    )

    candidates = stress[beg_idx + 1 :]
    if len(candidates) < nslices:
        raise RuntimeError(
            colors.FAIL
            + f"Only {len(candidates)} dump(s) after the transient cut, need {nslices}"
            + colors.ENDC
        )

    t_last = candidates[-1][1]
    targets = t_beg + (np.arange(1, nslices + 1) / nslices) * (t_last - t_beg)
    times = np.array([t for _, t, _ in candidates])

    slices: List[Slice] = []
    for k, target in enumerate(targets, start=1):
        folder, t_end, dumps_end = candidates[int(np.argmin(np.abs(times - target)))]
        slices.append(
            Slice(k, folder, t_beg, t_end, dumps_beg, dumps_end, dt, sample_frequency)
        )

    # Nearest-dump snapping can collide if the dumps are sparser than the
    # requested slicing; keep the windows strictly distinct.
    seen, unique = set(), []
    for s in slices:
        if s.folder in seen:
            print(colors.WARNING + f"Dropping duplicate window {s.folder}" + colors.ENDC)
            continue
        seen.add(s.folder)
        unique.append(s)
    return beg_folder, unique


# ---------------------------------------------------------------------------
# remote work
# ---------------------------------------------------------------------------
def run_combine_avg(
    dir_local: str,
    beg_folder: str,
    slices: List[Slice],
    mesh: str = "mesh.xml",
    host: str = "hpc",
    user: str = "vb824",
) -> None:
    """One SSH call producing every ``convAvg_kNN.fld`` window."""
    dir_remote = get_remote_dir(dir_local)
    cmds = []
    for s in slices:
        cmds.append(f"rm -rf {s.fld}")
        cmds.append(
            f"FieldConvert -m combineAvg:fromfld={beg_folder} "
            f"{mesh} {s.folder} {s.fld}"
        )

    print(
        colors.OKBLUE
        + f"Running {len(slices)} combineAvg window(s) on {host}..."
        + colors.ENDC
    )
    remote_cmd = f"""
    source /etc/profile
    source ~/.bashrc
    set -e
    cd {dir_remote}
    {'; '.join(cmds)}
    """
    run_cmd(f'ssh {user}@{host} "{remote_cmd}"')
    print(colors.OKGREEN + "combineAvg windows completed." + colors.ENDC)


def flds2datapts(
    dir_local: str,
    slices: List[Slice],
    pts_base: str,
    mesh: str = "mesh.xml",
    host: str = "hpc",
    user: str = "vb824",
) -> None:
    """Interpolate every window .fld onto the *same* points file in one SSH call.

    This is the transpose of ``pp.ssh.fieldconvert.fld2datapts_batch`` (one field,
    many point sets); here it is many fields, one point set. It also deliberately
    does not ``rm *.dat`` in the output directory the way ``fld2datapts`` does,
    since all the slices share that directory.
    """
    dir_remote = get_remote_dir(dir_local)
    out_dir_rel = os.path.dirname(slices[0].datbase)

    run_cmd(
        f"rsync -av --relative {dir_local}/./{pts_base}.pts "
        f"{user}@{host}:{dir_remote}/"
    )

    cmds = [f"mkdir -p {out_dir_rel}"]
    for s in slices:
        cmds.append(f"rm -f {s.datbase}.dat")
        cmds.append(
            f"FieldConvert -m interppoints:fromxml={mesh}:fromfld={s.fld}"
            f":topts={pts_base}.pts {s.datbase}.dat"
        )

    print(
        colors.OKBLUE
        + f"Interpolating {len(slices)} window(s) onto {pts_base}.pts..."
        + colors.ENDC
    )
    remote_cmd = f"""
    source /etc/profile
    source ~/.bashrc
    set -e
    cd {dir_remote}
    {'; '.join(cmds)}
    """
    run_cmd(f'ssh {user}@{host} "{remote_cmd}"')

    os.makedirs(os.path.join(dir_local, out_dir_rel), exist_ok=True)
    includes = " ".join([f"--include='{os.path.basename(s.datbase)}.dat'" for s in slices])
    run_cmd(
        f"rsync -av {includes} --exclude='*' "
        f"{user}@{host}:{dir_remote}/{out_dir_rel}/ {dir_local}/{out_dir_rel}/"
    )
    print(colors.OKGREEN + "Interpolation completed." + colors.ENDC)


def ensure_points_file(dir_local: str, case: str, pts_base: str) -> None:
    """Reuse the production .pts if it is already there; the slices must be
    sampled on exactly the same grid as the production curves for the
    comparison to mean anything."""
    pts_path = os.path.join(dir_local, pts_base + ".pts")
    if os.path.isfile(pts_path):
        print(colors.OKBLUE + f"Reusing existing points file: {pts_path}" + colors.ENDC)
        return
    d, w = extract_depth_width(case)
    gen_points_file(
        ymin=-d, ymax=150, n=NPOINTS_Y, width=w, dir_local=dir_local,
        output_file=pts_base,
    )


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------
def plot_convergence(
    slices: List[Slice],
    curves: List[Tuple[np.ndarray, np.ndarray]],
    dir_local: str,
    case: str,
    dt: float,
    sample_frequency: int,
    out_pdf: str,
) -> None:
    d, w = extract_depth_width(case)
    nsteps = np.array([s.nsteps for s in slices], dtype=float)

    # Linear colour map in the number of time steps: the visual question is
    # whether the curves pile up at the high-steps end of the map.
    norm = Normalize(vmin=nsteps.min(), vmax=nsteps.max())
    cmap = plt.get_cmap("viridis")
    sm = ScalarMappable(norm=norm, cmap=cmap)

    fig = plt.figure(figsize=(14, 6))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.1, 1], hspace=0.35, wspace=0.25)
    ax = fig.add_subplot(gs[:, 0])
    ax_win = fig.add_subplot(gs[0, 1])
    ax_dev = fig.add_subplot(gs[1, 1])

    x_ref, N_ref = curves[-1]

    for s, (x, Nx) in zip(slices, curves):
        ax.plot(x, Nx, color=cmap(norm(s.nsteps)), lw=1.2)
    # The longest window is the production curve: mark it so the eye has a target.
    ax.plot(x_ref, N_ref, color="k", lw=1.6, ls="--",
            label=f"longest window (T_avg = {slices[-1].T_avg:.0f})")

    ax.set_xlabel("x")
    ax.set_ylabel("N")
    ax.set_title(
        f"d = {d:g}, w = {w:g}  |  dt = {dt:g} "
        f"(sampled every {sample_frequency} steps)\n"
        f"transient removed: t > {slices[0].t_beg:.0f}"
    )
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.3)

    cbar = fig.colorbar(sm, ax=ax, pad=0.01)
    cbar.set_label("time steps used in the average")

    # --- clustering diagnostics --------------------------------------------
    # N averaged over the same downstream window used by computeDeltaN, i.e.
    # the single number that ends up in the DeltaN maps.
    x_start, x_end = w + 80.0, w + 130.0
    Nbar, dev = [], []
    for x, Nx in curves:
        idx = np.where((x >= x_start) & (x <= x_end))[0]
        Nbar.append(np.mean(Nx[idx]) if idx.size else np.nan)
        # Deviation from the longest window on their common x range.
        lo, hi = max(x[0], x_ref[0]), min(x[-1], x_ref[-1])
        xc = x_ref[(x_ref >= lo) & (x_ref <= hi)]
        dev.append(np.max(np.abs(np.interp(xc, x, Nx) - np.interp(xc, x_ref, N_ref))))
    Nbar, dev = np.array(Nbar), np.array(dev)

    ax_win.scatter(nsteps, Nbar, c=nsteps, cmap=cmap, norm=norm,
                   edgecolors="k", zorder=3, s=55)
    ax_win.axhline(Nbar[-1], color="k", ls="--", lw=1)
    ax_win.fill_between(
        [nsteps.min(), nsteps.max()], Nbar[-1] - 0.05, Nbar[-1] + 0.05,
        color="k", alpha=0.08, label="+/- 0.05 in N",
    )
    ax_win.set_ylabel(f"N averaged over\n{x_start:.0f} < x < {x_end:.0f}")
    ax_win.set_title("clustering of the DeltaN window value", fontsize=10)
    ax_win.legend(fontsize=8)
    ax_win.grid(alpha=0.3)

    ax_dev.scatter(nsteps[:-1], dev[:-1], c=nsteps[:-1], cmap=cmap, norm=norm,
                   edgecolors="k", zorder=3, s=55)
    ax_dev.set_yscale("log")
    ax_dev.set_xlabel("time steps used in the average")
    ax_dev.set_ylabel("max |N - N_longest|")
    ax_dev.grid(alpha=0.3, which="both")

    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
    print(colors.OKGREEN + f"Figure written to {out_pdf}" + colors.ENDC)

    print(colors.OKBLUE + "\nConvergence summary" + colors.ENDC)
    print(f"{'k':>3s} {'T_avg':>10s} {'steps':>12s} {'N_window':>10s} {'max|dN|':>10s}")
    for s, nb, dv in zip(slices, Nbar, dev):
        print(f"{s.k:3d} {s.T_avg:10.1f} {s.nsteps:12d} {nb:10.4f} {dv:10.4f}")


# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir", default=DEFAULT_DIR, help="case directory")
    parser.add_argument("-k", "--nslices", type=int, default=NSLICES)
    parser.add_argument("--mintime", type=float, default=MINTIME)
    parser.add_argument("--force", action="store_true",
                        help="redo the remote work even if the .dat files exist")
    parser.add_argument("--no-remote", action="store_true",
                        help="skip all remote work and only replot existing .dat files")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    dir_local = os.path.abspath(args.dir).rstrip("/")
    case = os.path.basename(dir_local)
    pts_base = f"data/pointsavg_n{NPOINTS_Y}"

    dt, sample_frequency = read_run_parameters(dir_local)
    dir_mounted = get_mounted_dir(dir_local)
    stress = scan_stress_folders(dir_mounted)
    print(colors.OKBLUE + f"Found {len(stress)} stress dump(s) in {case}" + colors.ENDC)

    beg_folder, slices = build_slices(
        stress, args.mintime, args.nslices, dt, sample_frequency
    )
    for s in slices:
        print("  " + repr(s))

    os.makedirs(os.path.join(dir_local, CONV_SUBDIR), exist_ok=True)
    missing = [s for s in slices
               if not os.path.isfile(os.path.join(dir_local, s.datbase + ".dat"))]

    if args.no_remote:
        if missing:
            raise SystemExit(
                colors.FAIL
                + f"--no-remote given but {len(missing)} .dat file(s) are missing"
                + colors.ENDC
            )
    elif args.force or missing:
        todo = slices if args.force else missing
        patch_start_folder(dir_mounted, beg_folder)
        run_combine_avg(dir_local, beg_folder, todo)
        ensure_points_file(dir_local, case, pts_base)
        flds2datapts(dir_local, todo, pts_base)
    else:
        print(colors.OKBLUE + "All slice .dat files present, skipping remote work"
              + colors.ENDC)

    with open(os.path.join(dir_local, CONV_SUBDIR, "slices.json"), "w") as f:
        json.dump(
            {"case": case, "dt": dt, "sample_frequency": sample_frequency,
             "start_folder": beg_folder, "mintime": args.mintime,
             "slices": [{"k": s.k, "folder": s.folder, "t_beg": s.t_beg,
                         "t_end": s.t_end, "T_avg": s.T_avg,
                         "nsamples": s.nsamples, "nsteps": s.nsteps,
                         "fld": s.fld, "dat": s.datbase + ".dat"} for s in slices]},
            f, indent=2,
        )

    curves, kept = [], []
    for s in slices:
        x, Nx = computeNx(os.path.join(dir_local, s.datbase + ".dat"), DO_LOO)
        if len(x) == 0:
            print(colors.WARNING + f"Skipping slice k={s.k}: no data" + colors.ENDC)
            continue
        curves.append((x, Nx))
        kept.append(s)
    if not curves:
        raise SystemExit(colors.FAIL + "No N factor curve could be computed" + colors.ENDC)

    out_pdf = os.path.join(dir_local, CONV_SUBDIR, f"{CONV_PREFIX}_nfactor_{case}.pdf")
    plot_convergence(kept, curves, dir_local, case, dt, sample_frequency, out_pdf)
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
