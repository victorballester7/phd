"""
N factor curves of the two gap families from the long linearised runs, with a
statistical error bar, against the short-window global curves used so far:

    d1.5    constant depth, w = 10 ... 50
    w40     constant width, d = 0.25 ... 1.75 (d1.5_w40 is shared)

It replaces the former ``plots/nfactorCurvesLocal.py``, which drew the same
figure from the short runs. There the global curve of a case was n(x) read off
one averaged field: the last cumulative stress file minus the
one at t ~ 5000 (``pp.ssh.fieldconvert.combineAvg``), i.e. a single window of
about 18000 time units (43000 for d1.75_w40). Here the same cases were run for
7.3e4 to 1.05e5 time units and the averaging is redone from every stress file
the run wrote.

How the average is built
------------------------
The ReynoldsStresses filter writes *cumulative* averages: mesh_k_stress.fld is
the average over [0, T_k], T_k = 1400 k, built from N_k samples. The amplitude
only needs

    Q = <u>^2 + <v>^2 + <u'u'> + <v'v'> = <u^2 + v^2>,

a raw second moment, which is linear in the samples. The average of Q over any
window (T_i, T_j] is therefore exact from two files,

    Q_ij = (N_j Q_j - N_i Q_i) / (N_j - N_i),

and so is S(x) = int Q dy = A(x)^2, because the y integral is linear too. This
is what combineAvg does with a negative NumberOfFieldDumps, without having to
rewrite any Info.xml on the cluster.

Every stress file was interpolated once onto the same points as the old
``data/pointsavg_n600.pts`` and reduced to Q(x, y) on the cluster
(``reduce_one.py`` in the run directory), giving data/nconv/Q_k.npz here.

The representative curve and its error bar
------------------------------------------
    n(x)        n of the whole post-transient window (T_START, T_end]: the best
                estimate the run can give.
    sigma_n(x)  standard error of n(x) by batch means. The window is cut into M
                consecutive blocks of BLOCK_FILES files each (BLOCK_FILES * 1400
                time units). Each block gives its own n_b(x), referenced to the
                same station x0 as n(x), and

                    sigma_n = std_b(n_b) / sqrt(M).

                The blocks are long against the correlation time of the forced
                signal, so the n_b are close to independent -- check the lag-1
                autocorrelation printed per case -- and, unlike the Isserlis
                estimate of ``computeNx(return_sigma=True)``, nothing has to be
                assumed about the integral time scale.

sigma_n is printed and saved, not drawn. The block-to-block spread std_b(n_b) scaled
to the old window length, std_b * sqrt(L_block / T_old), is printed as well:
it is how far a single window as long as the old one can be expected to fall
from the converged value.

The long d1.5 runs use a new mesh (``mesh.xml``; the old one is
``meshold.xml``) and every long run uses dt = 0.007 instead of 0.004 / 0.006,
so the old/new difference is not averaging time alone. (The w40 meshes are
regenerated from the same .geo as the old ones, so there only dt changed.) To
separate the two, the new runs are also averaged over a window as long as the
old one, (T_START, T_START + T_old] -- the n_same column of the printed table.

Figure, <fam> = d1.5 or w40, in images/ and copied to the paper:
    nfactorCurvesLocal<fam>.pdf
        global n (long average) solid, local LSA dashed, flat plate always and
        the backward facing step in d1.5. The local curve is the incoherent sum
        over the scanned frequency band (``pp.localNfactor.synthesise``, flat
        weights) from the cached Orr-Sommerfeld scans of
        ``scripts/predictNfactorLocal.py``. The flat plate and the step were not
        rerun; their global curves come from their averaged field.

Run from scripts/pp: ``uv run python plots/nfactorCurvesLongAvg.py [--sync]``;
``--sync`` first copies data/nconv/ from the cluster.
"""

import glob
import os
import shutil
import subprocess
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, NullLocator

from pp.DeltaN_computation import computeNx
from pp.fileManagement import extract_depth_width
from pp.figureFrame import FigureFrame
from pp.localNfactor import ModeSet, synthesise

plt.style.use("plots/style/jfm.mplstyle")

RE = 1000
BLASIUS_C = 1.7207876573

ROOT = "/home/victor/Desktop/PhD"
SRC_ROOT = os.path.join(ROOT, "src")
CASE_ROOT = os.path.join(SRC_ROOT, f"incGapRe{RE}", "directLinearSolver", "blowingSuction")
REMOTE_ROOT = f"vb824@hpc:Desktop/PhD/runs/incGapRe{RE}/directLinearSolver/blowingSuction"
NSAMPLE = 600
DATA_ROOT = os.path.join(ROOT, "data", "localNfactor")

# the cached local scan: F = 10 to 500 in 28 steps, 81 collocation nodes, free
# stream at y = 45, stations out to x = 700
TAG = "F10-500x28_n81_y45_X700"

COL_FLAT = "tab:blue"
COL_BFS = "tab:orange"

FAMILIES = {
    "d1.5": dict(
        cases={
            "d1.5_w10": r"$w=10$",
            "d1.5_w15": r"$w=15$",
            "d1.5_w20": r"$w=20$",
            "d1.5_w30": r"$w=30$",
            "d1.5_w40": r"$w=40$",
            "d1.5_w50": r"$w=50$",
        },
        cmap="Greens",
        extras=[("bfs_d1.5", "BFS", COL_BFS)],
    ),
    "w40": dict(
        cases={
            "d0.25_w40": r"$d=0.25$",
            "d0.5_w40": r"$d=0.50$",
            "d0.75_w40": r"$d=0.75$",
            "d1_w40": r"$d=1.00$",
            "d1.25_w40": r"$d=1.25$",
            "d1.5_w40": r"$d=1.50$",
            "d1.75_w40": r"$d=1.75$",
        },
        cmap="Reds",
        extras=[],
    ),
}
ALL_CASES = list(dict.fromkeys(n for fam in FAMILIES.values() for n in fam["cases"]))

# the short runs these are compared against were renamed *_normalTime when the
# long ones took their directory
OLD_SUFFIX = "_normalTime"

# averaging
T_START = 5000.0  # transient cut, as in the old combineAvg: first file with T >= this
T_OLD = 18000.0  # window of the old curves, (6000, 24000] ...
T_OLD_CASE = {"d1.75_w40": 42720.0}  # ... except this one, (5280, 48000]
BLOCK_FILES = 5  # 7000 time units per batch
YMAX = 150.0  # top of the L2 integral, as computeAmplitude

# window the n factor is averaged over, as computeDeltaN, ticked on the global
# curve of each gap case in a darkened shade of its colour
AVG_OFFSET = 80.0
AVG_LENGTH = 50.0
AVG_NMARK = 6
AVG_DARKEN = 0.85

# frame
YLIM = (0.0, 7.5)
YTICKS = [0, 2, 4, 6]
XLIM_X = (-80.0, 600.0)
RE_TICKS = [1000, 1200, 1400, 1600]
X_TICKS = [-70, 0, 40, 200, 400, 600]

LW_GLOBAL = 1.0
LW_LOCAL = 0.9
LS_GLOBAL = "-"
LS_LOCAL = "--"

# gap shading, drawn behind the curves
SHADE_GREY = "0.88"  # a single gap
# nested gaps: grey levels bracketing SHADE_GREY, from the stretch every gap
# covers (lightest) to the one only the widest gap reaches (darkest)
SHADE_LIGHT = 0.95
SHADE_DARK = 0.83

NSIG = 2.0  # error bars of the printed table


def x_to_re(x):
    return RE * np.sqrt(np.maximum(1.0 + BLASIUS_C**2 * np.asarray(x) / RE, 0.0))


def re_to_x(r):
    return ((np.asarray(r) / RE) ** 2 - 1.0) * RE / BLASIUS_C**2


def gap_width(name):
    """Streamwise extent of the gap, inf for a step and 0 for the flat plate."""
    if name == "flat":
        return 0.0
    if name.startswith("bfs_"):
        return np.inf  # a step has no downstream edge
    return float(extract_depth_width(name)[1])


def t_old(name):
    return T_OLD_CASE.get(name, T_OLD)


def data_dir(name):
    return os.path.join(CASE_ROOT, name, "data", "nconv")


def sync_data():
    for name in ALL_CASES:
        os.makedirs(data_dir(name), exist_ok=True)
        subprocess.run(
            f"rsync -a --include='Q_*.npz' --exclude='*' {REMOTE_ROOT}/{name}/nconv/ {data_dir(name)}/",
            shell=True,
            check=True,
        )


# ---------------------------------------------------------------------------
# the cumulative record of a case
# ---------------------------------------------------------------------------
class Record:
    """S_k(x) = int_0^YMAX Q_k dy of every cumulative stress file of a run."""

    def __init__(self, name):
        files = glob.glob(os.path.join(data_dir(name), "Q_*.npz"))
        if not files:
            raise FileNotFoundError(f"no Q_*.npz in {data_dir(name)}; run with --sync")
        files.sort(key=lambda f: int(os.path.basename(f)[2:-4]))
        T, N, S = [], [], []
        for f in files:
            z = np.load(f)
            if not T:
                self.x, y = z["x"], z["y"]
                # same truncation as computeAmplitude: points below index yindx
                jmax = np.argmin(np.abs(y - YMAX))
                self.y = y[:jmax]
            T.append(float(z["T"]))
            N.append(int(z["N"]))
            S.append(np.trapezoid(z["Q"][:, :jmax], self.y, axis=1))
        # file 0 is the empty average at t = 0
        self.T = np.concatenate(([0.0], T))
        self.N = np.concatenate(([0], N)).astype(float)
        self.S = np.vstack((np.zeros_like(S[0]), S))
        if np.any(np.diff(self.N) <= 0):
            raise ValueError(f"{name}: NumberOfFieldDumps not increasing")

    def index(self, t):
        """First file with T >= t."""
        return int(np.searchsorted(self.T, t - 1.0))

    def window(self, i, j):
        """S(x) averaged over (T_i, T_j]."""
        return (self.N[j] * self.S[j] - self.N[i] * self.S[i]) / (self.N[j] - self.N[i])


def ref_station(x, S):
    """Index of A0, the minimum amplitude upstream of the gap, as computeNx."""
    idx = np.where(x <= 0)[0]
    return int(np.argmin(S[idx]))


def nfactor(S, i0):
    """n = log(A / A0) with A = sqrt(S); the curve starts at the reference."""
    return 0.5 * np.log(S[i0:] / S[i0])


def window_mean(x, n, w):
    x0 = w + AVG_OFFSET
    return float(np.mean(np.interp(np.linspace(x0, x0 + AVG_LENGTH, AVG_NMARK), x, n)))


def analyse(name):
    rec = Record(name)
    w = gap_width(name)
    s, e = rec.index(T_START), len(rec.T) - 1

    S_full = rec.window(s, e)
    i0 = ref_station(rec.x, S_full)
    x = rec.x[i0:]
    n = nfactor(S_full, i0)

    # batch means, blocks laid back from the end so the last file is used
    M = (e - s) // BLOCK_FILES
    edges = e - BLOCK_FILES * np.arange(M, -1, -1)
    nb = np.array([nfactor(rec.window(a, b), i0) for a, b in zip(edges[:-1], edges[1:])])
    sd = np.std(nb, axis=0, ddof=1)
    se = sd / np.sqrt(M)

    # window value of every block, and its lag-1 autocorrelation as a check of
    # the independence the batch-means error bar rests on
    mb = np.array([window_mean(x, v, w) for v in nb])
    d = mb - mb.mean()
    rho1 = float(np.sum(d[1:] * d[:-1]) / np.sum(d * d))

    # the new run over the old window only, mesh effect without the extra time
    j_old = rec.index(rec.T[s] + t_old(name))
    n_same = nfactor(rec.window(s, j_old), i0)

    L_block = rec.T[edges[1]] - rec.T[edges[0]]
    return dict(
        name=name,
        w=w,
        x=x,
        n=n,
        se=se,
        sd_block=sd,
        M=M,
        L_block=L_block,
        T_start=rec.T[s],
        T_end=rec.T[e],
        mean=window_mean(x, n, w),
        mean_se=float(np.std(mb, ddof=1) / np.sqrt(M)),
        mean_sd_old=float(np.std(mb, ddof=1) * np.sqrt(L_block / t_old(name))),
        mean_same=window_mean(x, n_same, w),
        rho1=rho1,
    )


def old_curve(name):
    """
    Global n(x) read off one averaged field: the short run of a gap case, or
    the only run of the flat plate and of the step. None when there is none.
    """
    if name == "flat":
        p = os.path.join(SRC_ROOT, f"flatPlateRe{RE}inc", "directLinearSolver", "blowingSuction", "data", f"pointsavg_n{NSAMPLE}.dat")
    elif name.startswith("bfs_"):
        p = os.path.join(SRC_ROOT, f"bfsRe{RE}inc", "directLinearSolver", "blowingSuction", name[4:], "data", f"pointsavg_n{NSAMPLE}.dat")
    else:
        p = os.path.join(CASE_ROOT, name + OLD_SUFFIX, "data", f"pointsavg_n{NSAMPLE}.dat")
    if not os.path.isfile(p):
        return None
    x, n = computeNx(p, False)
    if not len(x):
        return None
    if name == "flat":
        # the last station of the flat plate run sits where the 2D
        # extrapolation of the averaged field breaks down; it is replaced by
        # the linear continuation of the two before it
        n[-1] = 2.0 * n[-2] - n[-3]
    return x, n


def local_curve(name):
    """n(x) predicted by local theory, or None when the case was not scanned."""
    path = os.path.join(DATA_ROOT, f"modes_Re{RE}_{name}_{TAG}.npz")
    if not os.path.isfile(path):
        return None
    z = np.load(path, allow_pickle=True)
    ms = ModeSet(
        x=z["x"],
        F=z["F"],
        omega=z["omega"],
        alpha=z["alpha"],
        n_global=int(z["n_global"]),
        meta=dict(z["meta"].item()),
    )
    return synthesise(ms, None)


def save_table(r):
    out = os.path.join(CASE_ROOT, r["name"], "data", "nfactor_longavg.dat")
    np.savetxt(
        out,
        np.column_stack((r["x"], r["n"], r["se"], r["sd_block"])),
        fmt="%.8e",
        header=(
            f"n factor of {r['name']}, averaged over ({r['T_start']:.0f}, {r['T_end']:.0f}]\n"
            f"se: batch-means standard error, {r['M']} blocks of {r['L_block']:.0f}\n"
            "sd_block: std of n over the blocks\n"
            "x n se sd_block"
        ),
    )


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def _ramp(fam):
    return plt.get_cmap(fam["cmap"])(np.linspace(0.3, 1.0, len(fam["cases"])))


def shade_gaps(ax, entries):
    """
    Mark where the gap is.

    One grey band when every case in the family has the same gap. When the
    widths differ they are nested, all starting at x = 0, and each case shades
    the stretch its own gap adds to the one before it, in a grey that darkens
    with the width: the narrowest gap lightest, the widest darkest. The
    shading then reads like the single band of the constant width family while
    every downstream edge stays visible as a step in the grey.

    A step has no downstream edge and so nothing to delimit; it is left out.
    """
    widths = sorted({w for _, _, _, w in entries if 0.0 < w < np.inf})
    if not widths:
        return
    if len(widths) == 1:
        ax.axvspan(*x_to_re(np.array([0.0, widths[0]])), color=SHADE_GREY, lw=0, zorder=0)
        return
    lo = 0.0
    for grey, w in zip(np.linspace(SHADE_LIGHT, SHADE_DARK, len(widths)), widths):
        ax.axvspan(*x_to_re(np.array([lo, w])), color=str(grey), lw=0, zorder=0)
        lo = w


def plot_global_local(tag, fam, results, save_path):
    """
    Global n solid (long average for the gap cases), local n dashed, one panel
    per family. The flat plate and the step were not rerun, so their global
    curve is read off their averaged field.
    """
    frame = FigureFrame(frame_w=5.3, aspect_ratio=0.75, pad_l=0.71, pad_b=2.7, pad_t=0.7)
    fig, ax = frame.fig, frame.ax
    ax.grid(False)
    ramp = _ramp(fam)
    long_avg = {r["name"]: r for r in results}

    # flat plate first and the step next, so the gap curves are drawn over them
    entries = [("flat", "Flat plate", COL_FLAT, gap_width("flat"))]
    entries += [(n, lbl, c, gap_width(n)) for n, lbl, c in fam["extras"]]
    entries += [(n, lbl, ramp[i], gap_width(n)) for i, (n, lbl) in enumerate(fam["cases"].items())]
    shade_gaps(ax, entries)

    for name, label, col, w in entries:
        loc = local_curve(name)
        if name in long_avg:
            r = long_avg[name]
            ax.plot(x_to_re(r["x"]), r["n"], ls=LS_GLOBAL, lw=LW_GLOBAL, color=col, label=label)
            xm = np.linspace(w + AVG_OFFSET, w + AVG_OFFSET + AVG_LENGTH, AVG_NMARK)
            ax.plot(x_to_re(xm), np.interp(xm, r["x"], r["n"]), ls="", marker="|", markersize=3,
                    color=np.asarray(col) * AVG_DARKEN)
            has_global = True
        else:
            g = old_curve(name)
            has_global = g is not None
            if has_global:
                ax.plot(x_to_re(g[0]), g[1], ls=LS_GLOBAL, lw=LW_GLOBAL, color=col, label=label)
        if loc is not None:
            # no label: the dashed line shares its colour with the solid one,
            # and the linestyle is keyed by the second legend instead
            ax.plot(x_to_re(loc[0]), loc[1], ls=LS_LOCAL, lw=LW_LOCAL, color=col,
                    label=None if has_global else label)
        if not has_global or loc is None:
            print(f"  {name}: missing {'global' if not has_global else 'local'} curve")

    ax.set_xticks(RE_TICKS)
    ax.set_yticks(YTICKS)
    ax.set_xlim(*x_to_re(np.array(XLIM_X)))
    ax.set_ylim(*YLIM)
    ax.set_xlabel(r"$\mbox{\textit{Re}}_\mathrm{loc}$", labelpad=2)
    ax.set_ylabel(r"$n$", rotation=0, labelpad=10)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="both", top=False, right=True)
    ax.tick_params(axis="y", which="minor", length=2)
    frame.axis_on_top()
    ax2 = ax.secondary_xaxis("top", functions=(re_to_x, x_to_re))
    ax2.set_xlabel(r"$x$", labelpad=3)
    ax2.set_xticks(X_TICKS)
    ax2.xaxis.set_minor_locator(NullLocator())

    style_legend = ax.legend(
        handles=[
            Line2D([], [], color="0.25", ls=LS_GLOBAL, lw=LW_GLOBAL, label="Global LSA"),
            Line2D([], [], color="0.25", ls=LS_LOCAL, lw=LW_LOCAL, label="Local LSA"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=2,
        frameon=False,
    )
    ax.add_artist(style_legend)
    handles, labels = ax.get_legend_handles_labels()
    order = sorted(range(len(labels)), key=lambda i: labels[i] == "BFS")  # the step last
    ax.legend([handles[i] for i in order], [labels[i] for i in order], loc="upper center",
              bbox_to_anchor=(0.5, -0.3), ncol=3, columnspacing=1.0, frameon=False)

    fig.savefig(save_path, format="pdf")
    plt.close(fig)
    print(f"Saved {os.path.normpath(save_path)}")


def main():
    if "--sync" in sys.argv:
        sync_data()
    script_path = os.path.dirname(os.path.abspath(__file__))
    out_dirs = [
        os.path.join(script_path, "../../../images"),
        os.path.join(script_path, "../../../latex/papers/jfm_incNS/Images"),
    ]

    for tag, fam in FAMILIES.items():
        results = [analyse(name) for name in fam["cases"]]

        print(f"\n{tag}: window-averaged n over x in [w+{AVG_OFFSET:.0f}, w+{AVG_OFFSET + AVG_LENGTH:.0f}]")
        print(f"{'case':10s} {'window':>16s} {'M':>3s} {'n_long':>7s} {'2se':>6s} {'rho1':>6s} {'n_same':>7s} {'2sd_old':>7s} {'n_old':>7s}")
        for r in results:
            save_table(r)
            o = old_curve(r["name"])
            n_old = window_mean(*o, r["w"]) if o is not None else np.nan
            print(
                f"{r['name']:10s} ({r['T_start']:6.0f},{r['T_end']:7.0f}] {r['M']:3d} {r['mean']:7.3f} {NSIG * r['mean_se']:6.3f} "
                f"{r['rho1']:6.2f} {r['mean_same']:7.3f} {NSIG * r['mean_sd_old']:7.3f} {n_old:7.3f}"
            )

        fname = f"nfactorCurvesLocal{tag}.pdf"
        first = os.path.join(out_dirs[0], fname)
        plot_global_local(tag, fam, results, first)
        for d in out_dirs[1:]:
            shutil.copyfile(first, os.path.join(d, fname))


if __name__ == "__main__":
    main()
