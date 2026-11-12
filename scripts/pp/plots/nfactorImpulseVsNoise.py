"""
N factor curves from one impulse-response run against the awgn noise runs.

The noise runs (``<case>``) force the linearised equations with the Gaussian vortex times
awgn and average the Reynolds stresses over t in [5040, 24000]. The impulse runs
(``scripts/makeImpulseCases.py``) start from the forcing shape with no forcing and
accumulate int_0^T q q^T dt, which is the stationary covariance for a forcing b xi(t)
with xi white in time:

    <case>_impulse100   vortex weight 100, the noise runs' envelope
    <case>_impulse      vortex weight 1 (resolved, radius 0.7)
    <case>_impulse100src  vortex weight 100 on the source-refined mesh (scripts/refineSourceMesh.py, 2026-11-10):
                        on the production meshes the w=100 vortex is sampled at a single quadrature column
                        (its initial field there has a relative L2 error of ~1); the refined patch resolves it

Both are interpolated onto the noise run's ``data/pointsavg_n600.pts``
(``scripts/impulseStatsToPoints.py``), and N(x), DeltaN follow ``computeNx`` /
``computeDeltaN`` unchanged. DeltaN is like-for-like: each kind of run is referred to the
flat plate run with the same forcing.

Figures: images/nfactorImpulseVsNoise.pdf (N(x) per case, thin: N(x) - N_flat(x) against the
flat plate with the same forcing) and
images/nfactorImpulseVsNoiseDiff.pdf (N_impulse - N_noise against x - w).
"""

import contextlib
import io
import os
import re

import matplotlib.pyplot as plt
import numpy as np

from pp.DeltaN_computation import computeDeltaN, computeNx

plt.style.use("plots/style/jfm.mplstyle")

N = 600
CASES = ["d0.5_w10", "d0.5_w90", "d1.25_w90", "d1.5_w25", "d2_w28", "d2.5_w15", "d4_w15"]
KINDS = [
    # suffix, label, colour, linestyle
    ("", "noise (awgn)", "black", "-"),
    ("_impulse100", r"impulse, $w_g=100$", "#2a78d6", "--"),
    ("_impulse", r"impulse, $w_g=1$", "#eb6834", "-."),
    ("_srcsum", r"$\sum$ impulse per node (noise, $T\to\infty$)", "#1baf7a", ":"),
    ("_impulse100src", r"impulse, $w_g=100$, refined source mesh", "#eda100", "-"),
]
DIFF_AX = {"_impulse100": 0, "_impulse100src": 1, "_impulse": 2}  # difference panel of each kind
CONV_TOL = 0.01  # N converged once it stays within this of the final curve


X_OUTLET = 1000  # the flat-plate pts reach the outlet, where the impulse runs jump by ~0.8 (boundary effect)


def nx(path):
    with contextlib.redirect_stdout(io.StringIO()):
        x, Nx = computeNx(path, False)
    keep = x < X_OUTLET
    return x[keep], Nx[keep]


def convergence_time(conv_file, x_max):
    """First dump time after which N(x <= x_max) stays within CONV_TOL of the final curve."""
    d = np.load(conv_file)
    t, x, A = d["t"], d["x"], d["A"]
    i0 = np.argmin(A[-1][x <= 0])
    Nt = np.log(A[:, i0:] / A[:, [i0]])
    err = np.abs(Nt - Nt[-1])[:, x[i0:] <= x_max].max(axis=1)
    bad = np.nonzero(err > CONV_TOL)[0]
    return t[bad[-1] + 1] if len(bad) else t[0], t[-1]


def wall_time(output_file):
    """Time-integration wall time [h] from Nektar's output.txt."""
    if not os.path.isfile(output_file):
        return np.nan
    m = re.findall(r"Time-integration\s*:\s*([0-9.]+)s", open(output_file).read())
    return float(m[-1]) / 3600 if m else np.nan


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    src = os.path.join(script_path, "../../../src")
    base = os.path.join(src, "incGapRe1000/directLinearSolver/blowingSuction")
    flat = os.path.join(src, "flatPlateRe1000inc/directLinearSolver/blowingSuction")

    flats = {}
    for suffix, *_ in KINDS:
        f = os.path.join(f"{flat}{suffix}", "data", f"pointsavg_n{N}.dat")
        flats[suffix] = nx(f) if os.path.isfile(f) else None

    fig, axs = plt.subplots(2, 4, figsize=(10, 5.4), sharex=False, sharey=False)
    axs = axs.ravel()
    figd, (axd1, axd1r, axd2, axp) = plt.subplots(1, 4, figsize=(15.5, 3.4), gridspec_kw={"width_ratios": [1, 1, 1, 0.75]})
    axd2.sharey(axd1)
    axd1r.sharey(axd1)
    daxs = (axd1, axd1r, axd2)
    rows = []
    for ax, case in zip(axs, CASES):
        d, w = map(float, re.match(r"d([0-9.]+)_w([0-9.]+)", case).groups())
        curves = {}
        for suffix, label, c, ls in KINDS:
            f = os.path.join(base, case + suffix, "data", f"pointsavg_n{N}.dat")
            if not os.path.isfile(f):
                continue
            x, Nx = nx(f)
            dN = computeDeltaN(w, x, Nx, *flats[suffix]) if flats[suffix] is not None else np.nan
            curves[suffix] = (x, Nx, dN)
            ax.plot(x, Nx, color=c, linestyle=ls, linewidth=1.2 if suffix else 1.8, label=label)
            if flats[suffix] is not None:  # N(x) - N_flat(x), flat plate with the same forcing
                xf, Nf = flats[suffix]
                m = (x >= xf[0]) & (x <= xf[-1])
                ax.plot(x[m], Nx[m] - np.interp(x[m], xf, Nf), color=c, linestyle=ls,
                        linewidth=0.7 if suffix else 1.0, alpha=0.6, label="_nolegend_")
        ax.axvspan(0, w, color="0.88", zorder=0)
        ax.axvspan(w + 80, w + 130, color="0.75", alpha=0.35, zorder=0)
        ax.set_title(rf"$d={d:g}$, $w={w:g}$", fontsize=9)
        ax.set_xlim(-70, 1000)

        if "" not in curves:
            continue
        xn, Nn, dNn = curves[""]
        row = [case, dNn]
        for suffix, _, c, ls in KINDS[1:]:
            if suffix not in curves:
                row += [np.nan] * 5
                continue
            x, Nx, dN = curves[suffix]
            xc = xn[(xn >= max(x[0], xn[0])) & (xn <= min(x[-1], xn[-1]))]
            diff = np.interp(xc, x, Nx) - np.interp(xc, xn, Nn)
            if suffix == "_srcsum":
                axd1.plot(xc - w, diff, color="black", linestyle="--", linewidth=1.4,
                          label=rf"{case}: $\sum$ per-node impulses")
            else:
                daxs[DIFF_AX[suffix]].plot(xc - w, diff, linewidth=1.0, label=case)
            win = (xc >= w + 80) & (xc <= w + 130)
            conv_file = os.path.join(base, case + suffix, "data", "convergence.npz")
            t_conv = convergence_time(conv_file, x[-1])[0] if os.path.isfile(conv_file) else np.nan
            out = os.path.join(base, case + suffix, "data", "output.txt")
            row += [dN, np.abs(diff[xc <= w + 130]).max(), np.interp(w + 500, xc, diff), t_conv, wall_time(out)]
        rows.append(row)

    # flat plate: last panel, and its impulse - noise difference with w = 0
    axf = axs[len(CASES)]
    for suffix, label, c, ls in KINDS:
        if flats[suffix] is not None:
            axf.plot(*flats[suffix], color=c, linestyle=ls, linewidth=1.2 if suffix else 1.8, label=label)
    axf.set_title("flat plate", fontsize=9)
    axf.set_xlim(-70, 1000)
    if flats[""] is not None:
        xn, Nn = flats[""]
        for suffix, *_ in KINDS[1:]:
            if flats[suffix] is None:
                continue
            x, Nx = flats[suffix]
            xc = xn[(xn >= max(x[0], xn[0])) & (xn <= min(x[-1], xn[-1]))]
            if suffix == "_srcsum":
                axd1.plot(xc, np.interp(xc, x, Nx) - np.interp(xc, xn, Nn), color="0.45", linestyle="--",
                          linewidth=1.4, label=r"flat plate: $\sum$ per-node impulses")
            else:
                daxs[DIFF_AX[suffix]].plot(
                    xc, np.interp(xc, x, Nx) - np.interp(xc, xn, Nn), color="0.45", linestyle=":",
                    linewidth=1.6, label="flat plate")
            print(f"flat plate {suffix}: N - N_noise at x = 100, 300, 600, 900:",
                  " ".join(f"{np.interp(v, x, Nx) - np.interp(v, xn, Nn):6.3f}" for v in (100, 300, 600, 900)))

    h, l = axs[3].get_legend_handles_labels()  # d1.5_w25 has every kind
    h.append(plt.Line2D([], [], color="0.3", linewidth=0.8, alpha=0.6))
    l.append(r"thin: $N - N_{\mathrm{flat}}$ (same forcing)")
    fig.legend(h, l, loc="lower center", ncol=len(h), frameon=False, fontsize=8)
    for a in axs[4:]:
        a.set_xlabel(r"$x/\delta^*_0$")
    for a in axs[::4]:
        a.set_ylabel(r"$N$")
    fig.tight_layout(rect=(0, 0.06, 1, 1))

    for a, title in ((axd1, r"impulse $w_g=100$ $-$ noise"), (axd1r, r"impulse $w_g=100$, refined source $-$ noise"),
                     (axd2, r"impulse $w_g=1$ $-$ noise")):
        a.axhline(0, color="black", linewidth=0.8)
        a.axvspan(80, 130, color="0.75", alpha=0.35, zorder=0)
        a.set_xlabel(r"$(x-w)/\delta^*_0$")
        a.set_title(title, fontsize=9)
    axd1.set_ylabel(r"$N_{\mathrm{impulse}} - N_{\mathrm{noise}}$")
    axd1.legend(fontsize=7, ncol=2, frameon=False)

    # DeltaN parity: each kind of run against the flat plate with the same forcing
    dn = np.array([r[1] for r in rows])
    cols = (2, 7, 12, 17)
    lim = [0, 1.05 * np.nanmax([dn.max()] + [np.nanmax([r[i] for r in rows]) for i in cols])]
    axp.plot(lim, lim, color="black", linewidth=0.8, zorder=0)
    for i, (suffix, label, c, _) in zip(cols, KINDS[1:]):
        axp.plot(dn, [r[i] for r in rows], {"_impulse100": "o", "_impulse": "s", "_impulse100src": "^"}.get(suffix, "D"), color=c,
                 markersize=5, markeredgecolor="white", markeredgewidth=0.8, linestyle="none", label=label)
    axp.set_xlim(lim)
    axp.set_ylim(lim)
    axp.set_aspect("equal")
    axp.set_xlabel(r"$\Delta n$, noise")
    axp.set_ylabel(r"$\Delta n$, impulse")
    axp.set_title(r"$\Delta n$ ($x\in[w+80, w+130]$)", fontsize=9)
    axp.legend(fontsize=7, frameon=False, loc="upper left")
    figd.tight_layout()

    out = os.path.join(script_path, "../../../images/nfactorImpulseVsNoise.pdf")
    fig.savefig(out)
    figd.savefig(out.replace(".pdf", "Diff.pdf"))
    print(f"Saved {out} and {out.replace('.pdf', 'Diff.pdf')}")

    hdr = "case        dN_noise | " + " | ".join(
        f"{k}: dN  max|dN(x)|<=w+130  dN(w+500)  t_conv  wall[h]" for k in ("w100", "w1", "node sum", "w100 refined"))
    print(hdr)
    for r in rows:
        print(f"{r[0]:11s} {r[1]:7.3f} | " + " | ".join(
            f"{r[i]:6.3f} {r[i+1]:9.3f} {r[i+2]:12.3f} {r[i+3]:8.0f} {r[i+4]:7.2f}" for i in cols))


if __name__ == "__main__":
    main()
