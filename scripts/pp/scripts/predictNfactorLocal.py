"""
Predict the N factor of a gap flow from local (parallel flow) stability theory.

At every streamwise station the base flow profile is frozen and treated as
parallel, the spatial Orr-Sommerfeld problem is solved for a band of real
frequencies, and the resulting alpha_i(x; omega) are integrated into one
N_j(x) per frequency. The perturbation is then modelled as the incoherent sum
of those modes, each carrying a unit L2 norm eigenfunction and a constant
weight a_j:

    A(x)^2 = sum_j a_j^2 exp(2 N_j(x)),    N(x) = 1/2 log(A^2 / A0^2)

See ``pp.localNfactor`` for why the modes add in energy rather than amplitude,
and why the unit norm convention makes the eigenfunctions drop out of the
answer altogether.

The prediction is compared against the DNS N factor of the same case, taken
from ``pp.DeltaN_computation.computeNx`` so that both curves use the same
amplitude definition and the same reference station.
"""

# Single threaded BLAS: every matrix here is a few hundred square, where the
# OpenBLAS thread pool costs far more than the factorisation it parallelises.
# The parallelism is taken across stations instead, inside track_modes.
import os

for _v in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_v, "1")

import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pp.colors import colors
from pp.DeltaN_computation import computeDeltaN, computeNx
from pp.fileManagement import extract_depth_width
from pp.localNfactor import (
    ModeSet,
    blasius_stations,
    fit_weights,
    jump_diagnostic,
    load_sections,
    synthesise,
    track_modes,
)

_SCRIPT_DIR = Path(__file__).resolve().parent
_ROOT = (_SCRIPT_DIR / ".." / ".." / "..").resolve()

# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------
REYNOLDS = [1000] 

# Only cases holding both the base flow sections (for the local analysis) and
# the averaged fields (for the DNS N factor) can be compared; the list is
# filtered against the disk at run time.
CASES = {
    800: ["d4_w12", "d0.5_w90"],
    1000: [
        # constant depth, varying width
        "d1.5_w10",
        "d1.5_w15",
        "d1.5_w20",
        "d1.5_w30",
        "d1.5_w40",
        "d1.5_w50",
        # constant width, varying depth (d1.5_w40 is shared with the family
        # above and is listed once)
        "d0.25_w40",
        "d0.5_w40",
        "d0.75_w40",
        "d1_w40",
        "d1.25_w40",
        "d1.75_w40",
    ],
    3000: ["d4_w12", "d0.5_w90"],
}

# Backward facing steps, named by depth alone. A step is the w -> infinity
# limit of the gap, so it is the curve the constant depth family walks towards
# as the width grows, and it is worth carrying through the same local analysis.
# It has no downstream edge, hence no gap window: it is excluded from the
# Delta N summary below for the same reason the flat plate is.
BFS_CASES = {
    1000: ["d1.5"],
}

# Reduced frequency band, F = omega * 1e6 / Re. F = 500 is roughly the top of
# what the gap shear layer puts out; Blasius branch I at these Reynolds
# numbers sits near F = 100, so the band covers both.
F_MIN, F_MAX, F_NUM = 10.0, 500.0, 28

N_CHEB = 81  # collocation nodes (see docs/localNfactor.md for the convergence)
YMAX = 45.0  # free stream truncation, in delta*_le
NSAMPLE = 600  # profile resolution tag of the section files
# Stations past this are dropped: the Delta N window ends at w + 130 and the
# plots at 600, so the tail of the domain only costs time.
XMAX_STATION = 700.0

# "flat" : every frequency weighted equally, the zeroth order model for the
#          white in time forcing the DNS applies
# "fit"  : non negative least squares weights against the DNS amplitude
WEIGHTS = "flat"
FIT_WINDOW = None  # (xmin, xmax) to restrict the fit, or None for all x

CACHE_DIR = _ROOT / "data" / "localNfactor"
IMAGE_DIR = _ROOT / "images" / "localNfactor"
OVERWRITE = False


@dataclass
class Case:
    re: int
    name: str  # "d1.5_w30", or "flat"
    d: float
    w: float
    section_dir: Path | None  # None -> analytical Blasius stations
    avg_file: Path | None  # None -> no DNS curve to compare against


def _base(re: int) -> Path:
    return _ROOT / "src" / f"incGapRe{re}" / "directLinearSolver" / "blowingSuction"


def _flat_base(re: int) -> Path:
    return (
        _ROOT / "src" / f"flatPlateRe{re}inc" / "directLinearSolver" / "blowingSuction"
    )


def _bfs_base(re: int) -> Path:
    return _ROOT / "src" / f"bfsRe{re}inc" / "directLinearSolver" / "blowingSuction"


def is_gap(case: "Case") -> bool:
    """True for the gap cases, i.e. the ones that have a downstream edge."""
    return case.name != "flat" and not case.name.startswith("bfs_")


def discover(re: int) -> list[Case]:
    """Build the case list, skipping anything whose inputs are not on disk."""
    flat_avg = _flat_base(re) / "data" / f"pointsavg_n{NSAMPLE}.dat"
    # the flat plate sections were never extracted, so the reference uses the
    # analytical similarity profile -- which is what the flat plate is anyway
    out = [Case(re, "flat", 0.0, 0.0, None, flat_avg if flat_avg.is_file() else None)]

    for name in CASES.get(re, []):
        sec = _base(re) / name / "data" / "xSections"
        avg = _base(re) / name / "data" / f"pointsavg_n{NSAMPLE}.dat"
        if not sec.is_dir() or not list(sec.glob(f"points_n{NSAMPLE}_x*.dat")):
            print(colors.WARNING + f"Re{re} {name}: no sections, skipped" + colors.ENDC)
            continue
        d, w = extract_depth_width(name)
        out.append(Case(re, name, d, w, sec, avg if avg.is_file() else None))

    # the step carries the "bfs_" prefix in its cache and case name so that it
    # cannot collide with a gap of the same depth
    for name in BFS_CASES.get(re, []):
        sec = _bfs_base(re) / name / "data" / "xSections"
        avg = _bfs_base(re) / name / "data" / f"pointsavg_n{NSAMPLE}.dat"
        if not sec.is_dir() or not list(sec.glob(f"points_n{NSAMPLE}_x*.dat")):
            print(
                colors.WARNING + f"Re{re} bfs {name}: no sections, skipped" + colors.ENDC
            )
            continue
        d, _ = extract_depth_width(name)
        out.append(Case(re, f"bfs_{name}", d, 0.0, sec, avg if avg.is_file() else None))
    return out


def _cache_path(case: Case) -> Path:
    """
    Cache file name. Every setting that changes the stored alpha(x, F) is in
    the name, so editing one of them writes a new file instead of silently
    reusing a result computed under different settings.
    """
    return CACHE_DIR / (
        f"modes_Re{case.re}_{case.name}_F{F_MIN:g}-{F_MAX:g}x{F_NUM}"
        f"_n{N_CHEB}_y{YMAX:g}_X{XMAX_STATION:g}.npz"
    )


def run_case(case: Case) -> ModeSet:
    """Scan the modes of one case, reusing the cache when it is there."""
    path = _cache_path(case)
    if path.is_file() and not OVERWRITE:
        z = np.load(path, allow_pickle=True)
        print(colors.OKBLUE + f"Re{case.re} {case.name}: cached" + colors.ENDC)
        return ModeSet(
            x=z["x"],
            F=z["F"],
            omega=z["omega"],
            alpha=z["alpha"],
            n_global=int(z["n_global"]),
            meta=dict(z["meta"].item()),
        )

    F = np.linspace(F_MIN, F_MAX, F_NUM)
    if case.section_dir is None:
        # same station spacing as the gap cases so the curves line up
        x = np.concatenate(
            (np.arange(-70.0, 0.0, 5.0), np.arange(0.0, XMAX_STATION + 1.0, 10.0))
        )
        stations = blasius_stations(x, case.re)
    else:
        stations = load_sections(
            str(case.section_dir), f"points_n{NSAMPLE}_x*.dat", xmax=XMAX_STATION
        )

    print(
        colors.OKBLUE
        + f"Re{case.re} {case.name}: {len(stations)} stations x {F_NUM} frequencies"
        + colors.ENDC
    )
    t0 = time.time()
    ms = track_modes(stations, float(case.re), F, n=N_CHEB, ymax=YMAX)
    print(
        colors.OKGREEN + f"Re{case.re} {case.name}: {time.time() - t0:.0f} s, "
        f"{np.mean(np.isfinite(ms.alpha)) * 100:.1f}% of (station, F) pairs "
        f"had an admissible mode, "
        f"{len(jump_diagnostic(ms))} branch jumps" + colors.ENDC
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        x=ms.x,
        F=ms.F,
        omega=ms.omega,
        alpha=ms.alpha,
        n_global=ms.n_global,
        meta=np.array(ms.meta, dtype=object),
    )
    return ms


def dns_curve(case: Case) -> tuple[np.ndarray, np.ndarray] | None:
    if case.avg_file is None:
        return None
    x, n = computeNx(str(case.avg_file), doLoo=False)
    return (x, n) if len(x) else None


# --------------------------------------------------------------------------
def main() -> None:
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    summary: dict[int, np.ndarray] = {}

    for re in REYNOLDS:
        cases = discover(re)
        if len(cases) < 2:
            print(colors.WARNING + f"Re{re}: nothing to do" + colors.ENDC)
            continue

        # one case at a time: track_modes already spreads its stations over
        # every core, so a second pool here would only oversubscribe
        modes = [run_case(c) for c in cases]

        flat_ms, flat_case = modes[0], cases[0]
        flat_dns = dns_curve(flat_case)
        x_flat_pred, n_flat_pred = synthesise(flat_ms, None)

        fig, (ax_p, ax_d) = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
        cmap = plt.get_cmap("tab10")
        rows: list[list[float]] = []

        for k, (case, ms) in enumerate(zip(cases, modes)):
            dns = dns_curve(case)
            weights = None
            if WEIGHTS == "fit" and dns is not None:
                weights = fit_weights(ms, dns[0], dns[1], xfit=FIT_WINDOW)

            xp, npred = synthesise(ms, weights)
            c = "black" if case.name == "flat" else cmap(k % 10)
            lw = 2.5 if case.name == "flat" else 1.5
            if case.name == "flat":
                lbl = "flat plate"
            elif case.name.startswith("bfs_"):
                lbl = f"BFS d={case.d:g}"
            else:
                lbl = f"d={case.d:g}, w={case.w:g}"

            ax_p.plot(xp, npred, color=c, lw=lw, label=lbl)
            if dns is not None:
                ax_d.plot(dns[0], dns[1], color=c, lw=lw, label=lbl)

            if not is_gap(case):
                continue

            # how much of the predicted Delta N is produced inside the gap
            xf, nfz = synthesise(ms, weights, freeze=(0.0, case.w))
            dn_pred = computeDeltaN(case.w, xp, npred, x_flat_pred, n_flat_pred)
            dn_froz = computeDeltaN(case.w, xf, nfz, x_flat_pred, n_flat_pred)
            dn_dns = np.nan
            if dns is not None and flat_dns is not None:
                dn_dns = computeDeltaN(case.w, dns[0], dns[1], flat_dns[0], flat_dns[1])
            rows.append([case.d, case.w, dn_pred, dn_froz, dn_dns])

            ax_p.plot(xf, nfz, color=c, lw=1, ls=":")

        for a, t in (
            (ax_p, "local theory (dotted: gap contribution removed)"),
            (ax_d, "DNS"),
        ):
            a.set_title(f"Re = {re}: {t}", fontsize=10)
            a.set_xlabel("x")
            a.grid(alpha=0.3)
            a.set_xlim(-80, 600)
        ax_p.set_ylabel("N(x)")
        ax_p.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(IMAGE_DIR / f"Nfactor_local_vs_dns_Re{re}.pdf")

        if rows:
            arr = np.array(rows)
            print(colors.HEADER + f"\n  Re = {re}:  Delta N" + colors.ENDC)
            print(
                f"  {'d':>6} {'w':>6} {'local':>9} {'no gap':>9} {'DNS':>9} {'err':>9}"
            )
            for d, w, a, b, c in arr:
                print(f"  {d:6.2f} {w:6.1f} {a:9.3f} {b:9.3f} {c:9.3f} {a - c:9.3f}")
            np.savetxt(
                CACHE_DIR / f"deltaN_local_vs_dns_Re{re}.dat",
                arr,
                header="d w deltaN_local deltaN_local_nogap deltaN_dns",
                comments="# ",
                fmt="%.6f",
            )
            summary[re] = arr

        # growth rate map and the per frequency fan, one figure per gap case
        for case, ms in zip(cases, modes):
            if not is_gap(case):
                continue
            fig2, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))

            g = ms.growth
            # Scale the colour map on the *amplified* side. The damped high F
            # modes reach -alpha_i ~ -0.2 and would otherwise flatten the whole
            # unstable region, which is the part worth seeing, to white.
            vm = float(np.nanmax(g))
            vm = max(vm, 1e-3)
            im = a1.pcolormesh(
                ms.x,
                ms.F,
                np.ma.masked_invalid(g),
                cmap="coolwarm",
                vmin=-vm,
                vmax=vm,
                shading="auto",
            )
            a1.contour(
                ms.x, ms.F, np.nan_to_num(g), levels=[0.0], colors="k", linewidths=1
            )
            fig2.colorbar(im, ax=a1, label=r"$-\alpha_i$", extend="min")
            a1.axvspan(0, case.w, color="grey", alpha=0.3)
            a1.set_xlabel("x")
            a1.set_ylabel("F")
            a1.set_xlim(-80, 600)
            a1.set_title(
                f"Re{re} d={case.d:g} w={case.w:g}: local growth rate", fontsize=10
            )

            N = ms.nfactor()
            for j in range(0, ms.F.size, max(1, ms.F.size // 10)):
                a2.plot(ms.x, N[j], lw=1, label=f"F={ms.F[j]:.0f}")
            xp, npred = synthesise(ms, None)
            a2.plot(xp, npred, "k-", lw=2.5, label="sum")
            a2.axvspan(0, case.w, color="grey", alpha=0.3)
            # the strongly damped modes run off to N ~ -80 and contribute
            # nothing to the sum; showing them would hide everything else
            a2.set_ylim(-8.0, max(2.0, float(np.nanmax(npred))) + 1.0)
            a2.set_xlabel("x")
            a2.set_ylabel(r"$N_j(x)$")
            a2.set_xlim(-80, 600)
            a2.grid(alpha=0.3)
            a2.legend(fontsize=7, ncol=2, loc="lower left")
            fig2.tight_layout()
            fig2.savefig(IMAGE_DIR / f"modes_Re{re}_{case.name}.pdf")
            plt.close(fig2)

    _summary_figure(summary)
    plt.show()


def _summary_figure(summary: dict[int, np.ndarray]) -> None:
    """
    Delta N from local theory against Delta N from the DNS, every case on one
    axis.

    The interesting number is not the scatter about the diagonal but the slope
    of the line through the origin: a single scalar, the same at both Reynolds
    numbers, means the local model has the physics and is simply missing a
    constant efficiency factor -- most plausibly the fraction of the gap shear
    layer's local growth that a disturbance of finite streamwise extent can
    actually collect while it crosses.
    """
    usable = {re: a for re, a in summary.items() if np.isfinite(a[:, 4]).any()}
    if not usable:
        return

    fig, ax = plt.subplots(figsize=(6, 6))
    allx, ally = [], []
    for marker, (re, a) in zip("os^v", sorted(usable.items())):
        m = np.isfinite(a[:, 4])
        ax.plot(a[m, 4], a[m, 2], marker, label=f"Re = {re}", markersize=7)
        allx.append(a[m, 4])
        ally.append(a[m, 2])

    x = np.concatenate(allx)
    y = np.concatenate(ally)
    slope = float(np.sum(x * y) / np.sum(x**2))
    resid = float(np.std(y - slope * x))

    lim = [0.0, 1.1 * max(x.max(), y.max())]
    ax.plot(lim, lim, "k--", lw=1, label="1:1")
    ax.plot(
        lim,
        [slope * v for v in lim],
        "r-",
        lw=1.5,
        label=f"slope {slope:.3f} (resid {resid:.2f})",
    )
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(r"$\Delta N$  (DNS)")
    ax.set_ylabel(r"$\Delta N$  (local theory, flat weights)")
    ax.set_title("Local stability prediction vs DNS")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(IMAGE_DIR / "deltaN_local_vs_dns.pdf")

    print(colors.HEADER + "\n  All cases:" + colors.ENDC)
    print(
        f"    slope through origin = {slope:.4f}  "
        f"(local theory over-predicts by {100 * (slope - 1):.1f}%)"
    )
    print(f"    residual scatter     = {resid:.3f} in Delta N")
    print(f"    correlation          = {np.corrcoef(x, y)[0, 1]:.5f}")


if __name__ == "__main__":
    main()
