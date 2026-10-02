"""
N factor curves of the d1.5 gap family from the long linearised runs, with a
statistical error bar, against the short-window global curves used so far.

Companion of ``plots/nfactorCurvesLocal.py``. There the global curve of a case
is n(x) read off one averaged field: the last cumulative stress file minus the
one at t ~ 5000 (``pp.ssh.fieldconvert.combineAvg``), i.e. a single window of
about 18000 time units. Here the same cases were run for 7.3e4 to 9.5e4 time
units and the averaging is redone from every stress file the run wrote.

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

The band drawn is n +- 2 sigma_n. The block-to-block spread std_b(n_b) scaled
to the old window length, std_b * sqrt(L_block / T_OLD), is printed as well:
it is how far a single 18000 unit window can be expected to fall from the
converged value.

The long runs use a new mesh (``mesh.xml``; the old one is ``meshold.xml``) and
dt = 0.007 instead of 0.006, so the old/new difference is not averaging time
alone. To separate the two, the new runs are also averaged over the old window
(T_START, T_START + T_OLD] -- the dotted curve of the convergence figure.

Figures (images/):
    nfactorCurvesLongAvgd1.5.pdf    n(x): long average +- 2 sigma vs the old curve
    nfactorConvergenced1.5.pdf      window-averaged n against averaging time
    nfactorWindowd1.5.pdf           window-averaged n against w, old vs new

Run from scripts/pp: ``uv run python plots/nfactorCurvesLongAvg.py [--sync]``;
``--sync`` first copies data/nconv/ from the cluster.
"""

import glob
import os
import subprocess
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import AutoMinorLocator, NullLocator

from pp.DeltaN_computation import computeNx
from pp.fileManagement import extract_depth_width
from pp.figureFrame import FigureFrame

plt.style.use("plots/style/jfm.mplstyle")

RE = 1000
BLASIUS_C = 1.7207876573

ROOT = "/home/victor/Desktop/PhD"
SRC_ROOT = os.path.join(ROOT, "src")
CASE_ROOT = os.path.join(SRC_ROOT, f"incGapRe{RE}", "directLinearSolver", "blowingSuction")
REMOTE_ROOT = f"vb824@hpc:Desktop/PhD/runs/incGapRe{RE}/directLinearSolver/blowingSuction"
NSAMPLE = 600

CASES = {
    "d1.5_w10": r"$w=10$",
    "d1.5_w15": r"$w=15$",
    "d1.5_w20": r"$w=20$",
    "d1.5_w30": r"$w=30$",
    "d1.5_w40": r"$w=40$",
    "d1.5_w50": r"$w=50$",
}
CMAP = "Greens"
COL_FLAT = "tab:blue"

# the short runs these are compared against were renamed *_normalTime when the
# long ones took their directory
OLD_SUFFIX = "_normalTime"

# averaging
T_START = 5000.0  # transient cut, as in the old combineAvg: first file with T >= this
T_OLD = 18000.0  # window of the old curves, (6000, 24000]
BLOCK_FILES = 5  # 7000 time units per batch
YMAX = 150.0  # top of the L2 integral, as computeAmplitude

# window the n factor is averaged over, as computeDeltaN / nfactorCurvesLocal.py
AVG_OFFSET = 80.0
AVG_LENGTH = 50.0
AVG_NMARK = 6

# frame, as nfactorCurvesLocal.py
YLIM = (0.0, 7.5)
YTICKS = [0, 2, 4, 6]
XLIM_X = (-80.0, 600.0)
RE_TICKS = [1000, 1200, 1400, 1600]
X_TICKS = [-70, 0, 40, 200, 400, 600]

LW = 1.0
LS_NEW = "-"
LS_OLD = "--"
BAND_ALPHA = 0.3
NSIG = 2.0
SHADE_ALPHA = 0.3


def x_to_re(x):
    return RE * np.sqrt(np.maximum(1.0 + BLASIUS_C**2 * np.asarray(x) / RE, 0.0))


def re_to_x(r):
    return ((np.asarray(r) / RE) ** 2 - 1.0) * RE / BLASIUS_C**2


def gap_width(name):
    return float(extract_depth_width(name)[1])


def data_dir(name):
    return os.path.join(CASE_ROOT, name, "data", "nconv")


def sync_data():
    for name in CASES:
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

    # running estimate: window value of (T_s, T_j] for every j
    run_T = rec.T[s + 1 : e + 1]
    run_m = np.array([window_mean(x, nfactor(rec.window(s, j), i0), w) for j in range(s + 1, e + 1)])

    # the new run over the old window only, mesh effect without the extra time
    j_old = rec.index(rec.T[s] + T_OLD)
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
        mean_sd_old=float(np.std(mb, ddof=1) * np.sqrt(L_block / T_OLD)),
        mean_same=window_mean(x, n_same, w),
        rho1=rho1,
        run_T=run_T,
        run_m=run_m,
        x_same=x,
        n_same=n_same,
        T_same=(rec.T[s], rec.T[j_old]),
    )


def old_curve(name):
    """The global curve as nfactorCurvesLocal.py drew it, from the short run."""
    if name == "flat":
        p = os.path.join(SRC_ROOT, f"flatPlateRe{RE}inc", "directLinearSolver", "blowingSuction", "data", f"pointsavg_n{NSAMPLE}.dat")
    else:
        p = os.path.join(CASE_ROOT, name + OLD_SUFFIX, "data", f"pointsavg_n{NSAMPLE}.dat")
    if not os.path.isfile(p):
        return None
    x, n = computeNx(p, False)
    if name == "flat":
        n[-1] = 2.0 * n[-2] - n[-3]  # as nfactorCurvesLocal.py
    return x, n


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
def _ramp():
    return plt.get_cmap(CMAP)(np.linspace(0.3, 1.0, len(CASES)))


def plot_curves(results, save_path):
    frame = FigureFrame(frame_w=5.3, aspect_ratio=0.75, pad_l=0.71, pad_b=3.0, pad_t=0.7)
    fig, ax = frame.fig, frame.ax
    ax.grid(False)
    ramp = _ramp()

    # gap staircase, as nfactorCurvesLocal.py
    lo = 0.0
    for col, r in zip(ramp, results):
        ax.axvspan(*x_to_re(np.array([lo, r["w"]])), color=col, alpha=SHADE_ALPHA, lw=0, zorder=0)
        lo = r["w"]

    flat = old_curve("flat")
    if flat is not None:
        ax.plot(x_to_re(flat[0]), flat[1], ls=LS_OLD, lw=LW, color=COL_FLAT, label="Flat plate")

    for col, r in zip(ramp, results):
        R = x_to_re(r["x"])
        ax.fill_between(R, r["n"] - NSIG * r["se"], r["n"] + NSIG * r["se"], color=col, alpha=BAND_ALPHA, lw=0)
        ax.plot(R, r["n"], ls=LS_NEW, lw=LW, color=col, label=CASES[r["name"]])
        old = old_curve(r["name"])
        if old is not None:
            ax.plot(x_to_re(old[0]), old[1], ls=LS_OLD, lw=0.8, color=col)
        x0 = r["w"] + AVG_OFFSET
        xm = np.linspace(x0, x0 + AVG_LENGTH, AVG_NMARK)
        ax.plot(x_to_re(xm), np.interp(xm, r["x"], r["n"]), ls="", marker="|", markersize=3, color=np.asarray(col) * 0.85)

    ax.set_xticks(RE_TICKS)
    ax.set_yticks(YTICKS)
    ax.set_xlim(*x_to_re(np.array(XLIM_X)))
    ax.set_ylim(*YLIM)
    ax.set_xlabel(r"$\mbox{\textit{Re}}_\ell$", labelpad=2)
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
            Line2D([], [], color="0.25", ls=LS_OLD, lw=0.8, label=r"$T_{\rm avg}=1.8\times10^4$"),
            (Patch(color="0.25", alpha=BAND_ALPHA, lw=0), Line2D([], [], color="0.25", ls=LS_NEW, lw=LW)),
        ],
        labels=[r"$T_{\rm avg}=1.8\times10^4$", rf"long average $\pm{NSIG:.0f}\sigma$"],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=2,
        frameon=False,
    )
    ax.add_artist(style_legend)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=3, columnspacing=1.0, frameon=False)

    fig.savefig(save_path, format="pdf")
    plt.close(fig)
    print(f"Saved {os.path.normpath(save_path)}")


def plot_convergence(results, save_path):
    """
    Window-averaged n of (T_start, T_start + T_avg] against T_avg, as a
    deviation from the long average: has it settled, and where did the old
    run land?
    """
    frame = FigureFrame(frame_w=5.3, aspect_ratio=0.75, pad_l=1.35, pad_r=0.3, pad_b=2.6, pad_t=0.2)
    fig, ax = frame.fig, frame.ax
    ax.grid(False)
    ax.axhline(0.0, color="0.5", lw=0.5)
    ax.axvline(T_OLD / 1e4, color="0.5", lw=0.5, ls=":")
    for col, r in zip(_ramp(), results):
        T = (r["run_T"] - r["T_start"]) / 1e4
        ax.plot(T, r["run_m"] - r["mean"], color=col, lw=LW, label=CASES[r["name"]])
        old = old_curve(r["name"])
        if old is not None:
            ax.plot([T_OLD / 1e4], [window_mean(*old, r["w"]) - r["mean"]], ls="", marker="o", mfc="none", color=col, ms=3.5)
    ax.set_xlim(0, None)
    ax.set_ylim(-0.12, 0.12)
    ax.set_xlabel(r"averaging time $T_{\rm avg}\,/10^4$", labelpad=2)
    ax.set_ylabel(r"$\bar n - \bar n_\infty$", labelpad=2)
    frame.jfm_ticks()
    frame.axis_on_top()
    handles = [Line2D([], [], color=c, lw=LW, label=CASES[r["name"]]) for c, r in zip(_ramp(), results)]
    handles += [
        Line2D([], [], ls="", marker="o", mfc="none", color="0.25", ms=3.5, label="old run"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, columnspacing=0.8, handlelength=1.2, handletextpad=0.4, frameon=False)
    fig.savefig(save_path, format="pdf")
    plt.close(fig)
    print(f"Saved {os.path.normpath(save_path)}")


def plot_window(results, save_path):
    """
    Window-averaged n against w, as a deviation from the long average:
    old run, new run over the old window, long average with its error bar.
    """
    frame = FigureFrame(frame_w=5.3, aspect_ratio=0.6, pad_l=1.35, pad_r=0.3, pad_b=2.6, pad_t=0.2)
    fig, ax = frame.fig, frame.ax
    ax.grid(False)
    w = np.array([r["w"] for r in results])
    m = np.array([r["mean"] for r in results])
    se = np.array([r["mean_se"] for r in results])
    same = np.array([r["mean_same"] for r in results])
    old = np.array([window_mean(*o, r["w"]) if (o := old_curve(r["name"])) is not None else np.nan for r in results])
    sd_old = np.array([r["mean_sd_old"] for r in results])

    ax.axhline(0.0, color="0.5", lw=0.5)
    ax.errorbar(w, 0 * m, yerr=NSIG * se, color="tab:green", marker="o", ms=3, lw=0, elinewidth=LW, capsize=2,
                label=rf"long average $\bar n_\infty \pm{NSIG:.0f}\sigma$")
    ax.errorbar(w + 0.8, same - m, yerr=NSIG * sd_old, color="0.45", marker="s", ms=2.5, lw=0, elinewidth=0.7, capsize=2,
                label=r"new run, $T_{\rm avg}=1.8\times10^4$, $\pm2\sigma_{18k}$")
    ax.plot(w - 0.8, old - m, ls="", marker="^", ms=3, color="tab:red", label=r"old run, $T_{\rm avg}=1.8\times10^4$")
    ax.set_xticks(w)
    ax.set_xlabel(r"$w$", labelpad=2)
    ax.set_ylabel(r"$\bar n - \bar n_\infty$", labelpad=2)
    frame.jfm_ticks()
    ax.xaxis.set_minor_locator(NullLocator())
    frame.axis_on_top()
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=1, frameon=False)
    fig.savefig(save_path, format="pdf")
    plt.close(fig)
    print(f"Saved {os.path.normpath(save_path)}")


def main():
    if "--sync" in sys.argv:
        sync_data()
    results = [analyse(name) for name in CASES]

    print(f"\nwindow-averaged n over x in [w+{AVG_OFFSET:.0f}, w+{AVG_OFFSET + AVG_LENGTH:.0f}]")
    print(f"{'case':10s} {'window':>15s} {'M':>3s} {'n_long':>7s} {'2se':>6s} {'rho1':>6s} {'n_same':>7s} {'2sd18k':>7s} {'n_old':>7s}")
    for r in results:
        save_table(r)
        o = old_curve(r["name"])
        n_old = window_mean(*o, r["w"]) if o is not None else np.nan
        print(
            f"{r['name']:10s} ({r['T_start']:6.0f},{r['T_end']:6.0f}] {r['M']:3d} {r['mean']:7.3f} {NSIG * r['mean_se']:6.3f} "
            f"{r['rho1']:6.2f} {r['mean_same']:7.3f} {NSIG * r['mean_sd_old']:7.3f} {n_old:7.3f}"
        )

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../images")
    plot_curves(results, os.path.join(out, "nfactorCurvesLongAvgd1.5.pdf"))
    plot_convergence(results, os.path.join(out, "nfactorConvergenced1.5.pdf"))
    plot_window(results, os.path.join(out, "nfactorWindowd1.5.pdf"))


if __name__ == "__main__":
    main()
