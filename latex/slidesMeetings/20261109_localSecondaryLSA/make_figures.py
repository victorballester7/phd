"""Figures for the 9 Nov 2026 meeting (local secondary stability of the d = 1.5 waves).

Run from anywhere: scripts/pp/.venv/bin/python make_figures.py
Uses data/stability/secondaryStabilityAll_d1.5.npz (scripts/pp/scripts/secondaryStabilityAll.py)
and the DNS profiles through scripts/pp/plots/secondaryVsDNS.py. Writes Images/*.png.
"""
import contextlib
import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PP = os.path.join(HERE, "../../../scripts/pp")
os.chdir(PP)  # the jfm style path in L2POvsX is relative to scripts/pp
sys.path.insert(0, "plots")
sys.path.insert(0, "scripts")

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

with contextlib.redirect_stdout(io.StringIO()):
    import secondaryVsDNS as S  # noqa: E402
from secondaryStabilityScan import load_wave  # noqa: E402

IMG = os.path.join(HERE, "Images")
os.makedirs(IMG, exist_ok=True)
D, ST, COL = S.D, S.ST, S.COL
SLAB = {"1e-3": r"$A=10^{-3}$", "1e-5": r"$A=10^{-5}$", "1e-9": r"$A=10^{-9}$", "0": r"$A=0$"}
SCOL = {"1e-5": "C1", "1e-9": "C2", "0": "C0"}


def save(fig, name):
    fig.savefig(f"{IMG}/{name}", dpi=250, bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


# ---- 1. the wave model at one station -----------------------------------------------------------
def fig_wave(w=65, x0=150.0):
    L = load_wave(w)
    W = L.at(x0)
    al = W["alpha"]
    m = np.abs(L.x - x0) <= 2 * np.pi / al
    X = L.x[m] - x0
    jy = int(np.argmax(np.abs(W["u"][1])))
    rec = lambda nmax: W["u"][0][jy].real + sum(2 * (W["u"][n][jy] * np.exp(1j * n * al * X)).real
                                                for n in range(1, nmax + 1))
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.2), gridspec_kw=dict(width_ratios=[1.3, 0.8, 1.0]))
    ax[0].plot(X, L.u[m, jy], "ko", ms=2.2, label="DNS snapshot")
    ax[0].plot(X, rec(1), "C3", label=r"$u_0+u_1$")
    ax[0].plot(X, rec(3), "C0--", label=r"$u_0+\dots+u_3$")
    ax[0].set_xlabel(r"$x-x_0$"); ax[0].set_ylabel(rf"$u'$ at $y={L.y[jy]:.1f}$")
    ax[0].legend(fontsize=6.5, loc="lower left", ncol=3, columnspacing=0.8, handlelength=1.5)
    ax[0].set_title(rf"(a) fit over $x_0\pm\lambda$, $w={w}$, $x_0-w={x0:g}$", fontsize=8)
    jj = L.y <= 12
    ax[1].plot(W["u"][0][jj].real, L.y[jj], "C2", label=r"$u_0$")
    for n, c in ((1, "C3"), (2, "C1"), (3, "C4")):
        ax[1].plot(2 * np.abs(W["u"][n][jj]), L.y[jj], c, label=rf"$2|u_{n}|$")
    ax[1].set_xlabel("velocity"); ax[1].set_ylabel(r"$y$"); ax[1].legend(fontsize=6.5)
    ax[1].set_title(rf"(b) profiles, $a={W['amp']:.3f}$", fontsize=8)
    for ww, x00, c in ((65, 150.0, COL[65]), (70, 150.0, COL[70]), (70, 400.0, "k")):
        Lw = L if ww == w else load_wave(ww)
        af = Lw.alpha_at(x00)
        ah = Lw.alpha_at(x00, method="hilbert")
        r = np.linspace(0.8, 1.2, 81)
        J = [np.mean(list(Lw._fit(x00, af * k, 3, 2.0)[1].values())) for k in r]
        ax[2].plot(r, J, color=c, label=rf"$w={ww}$, $x_0-w={x00:g}$")
        ax[2].plot(ah / af, np.mean(list(Lw._fit(x00, ah, 3, 2.0)[1].values())), "o", color=c, mfc="w")
    ax[2].set_xlabel(r"$\alpha/\alpha_{\rm best}$"); ax[2].set_ylabel("rel.\ rms residual")
    ax[2].set_ylim(0, 0.85); ax[2].legend(fontsize=6.5, loc="upper center", ncol=1)
    ax[2].set_title(r"(c) choosing $\alpha$ (open: phase-slope estimate)", fontsize=8)
    fig.tight_layout(w_pad=0.6)
    save(fig, "wave.png")


# ---- 2. all widths -----------------------------------------------------------------------------
def fig_all():
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.3))
    for w in S.W:
        iw = S.W.index(w)
        alive = D["alive"][iw]
        ax[0].plot(ST, np.where(alive, D["amp"][iw], np.nan), color=COL[w], label=rf"$w={w}$")
        ax[1].plot(ST, S.lsa(w), color=COL[w])
        xi, N = S.integrate(ST, np.maximum(np.nan_to_num(S.lsa(w)), 0), 0, 600)
        ax[2].plot(xi, N, color=COL[w])
    ax[0].set_yscale("log"); ax[0].set_ylabel(r"$a=\max_y 2|u_1|$"); ax[0].set_ylim(3e-3, 0.5)
    ax[0].legend(fontsize=6.5, ncol=2, loc="lower right")
    ax[0].set_title("(a) amplitude of the 2D wave", fontsize=8)
    ax[1].axhline(0, color="0.4", lw=0.6); ax[1].set_ylabel(r"$\sigma=\Re\lambda/c$")
    ax[1].set_title(r"(b) spatial growth rate, $\gamma=0$", fontsize=8)
    for A, lab in ((1e-17, r"$A\approx10^{-17}$ (no seed)"), (1e-9, r"$A=10^{-9}$"), (1e-5, r"$A=10^{-5}$")):
        v = np.log(1e-2 / A)
        ax[2].axhline(v, color="0.4", lw=0.6, ls=":")
        ax[2].text(12, v + 0.6, lab, ha="left", va="bottom", fontsize=6.5)
    ax[2].set_ylabel(r"$N=\int\sigma\,{\rm d}x$")
    ax[2].set_title(r"(c) $N$ from the gap; lines: $\ln(10^{-2}/A)$", fontsize=8)
    for a in ax:
        a.set_xlabel(r"$x-w$"); a.set_xlim(0, 600)
    fig.tight_layout(w_pad=0.6)
    save(fig, "all.png")


# ---- 3. against the DNS ------------------------------------------------------------------------
def fig_dns():
    """ln||zeta|| of the DNS and the integrated LSA growth (the sigma panels were dropped 2026-11-09)."""
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.4), sharey=True)
    for j, w in enumerate((65, 66, 70)):
        sb = S.lsa(w)
        for seed in ("1e-5", "1e-9", "0"):
            x, m, x_lo, x_hi, xs, s, e, n = S.dns(w, seed)
            c = SCOL[seed]
            ax[j].plot(x, m, color=c, lw=0.9)
            xi, N = S.integrate(ST, sb, x_lo, x_hi + 40)
            ax[j].plot(xi, np.interp(x_lo, x, m) + N, color=c, ls="--", lw=1.0)
        ax[j].axhline(np.log(1e-2), color="0.4", lw=0.6, ls=":")
        ax[j].set_ylim(-37, 1); ax[j].set_xlim(0, 500)
        ax[j].set_title(rf"$w={w}$", fontsize=8.5); ax[j].set_xlabel(r"$x-w$")
    ax[0].set_ylabel(r"$\ln\|\zeta\|$")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=SCOL[k], label=f"DNS {SLAB[k]}") for k in ("1e-5", "1e-9", "0")]
    h += [Line2D([], [], color="0.4", ls="--", label=r"DNS$(x_{\rm lo})+\int\sigma_{\rm LSA}{\rm d}x$")]
    fig.legend(handles=h, loc="upper center", ncol=4, fontsize=6.5, bbox_to_anchor=(0.5, 1.07))
    fig.tight_layout(w_pad=0.4)
    save(fig, "dns.png")


if __name__ == "__main__" and len(sys.argv) == 1:
    fig_wave()
    fig_all()
    fig_dns()


# ---- 4. growth rates only (meeting version) -------------------------------------------------------
def fig_sigma():
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.5), sharey=True)
    rows = []
    for j, w in enumerate((65, 66, 70)):
        sb, sa = S.lsa(w), S.lsa(w, box=False)
        ax[j].plot(ST, sb, "k-", lw=1.3)
        ax[j].plot(ST, sa, "k:", lw=1.0)
        d = {}
        for seed in ("1e-5", "1e-9", "0"):
            x, m, x_lo, x_hi, xs, s, e, n = S.dns(w, seed)
            d[seed] = (xs, s)
            ax[j].errorbar(xs, s, e, fmt="o", color=SCOL[seed], ms=2.5, capsize=1.5, elinewidth=0.6)
        for X in (100, 150, 200, 300, 400):
            g = lambda k: (np.interp(X, *d[k]) if d[k][0].size and d[k][0][0] - 1 <= X <= d[k][0][-1] + 1 else np.nan)
            rows.append((w, X, np.interp(X, ST, sb), np.interp(X, ST, sa), g("1e-9"), g("0")))
        ax[j].set_title(rf"$w={w}$", fontsize=8.5); ax[j].set_xlabel(r"$x-w$")
        ax[j].set_xlim(0, 500); ax[j].set_ylim(0, 0.16)
    ax[0].set_ylabel(r"spatial growth rate $\sigma$")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color="k", lw=1.3, label=r"LSA, $\beta$ of the $L_z=6$ box"),
         Line2D([], [], color="k", ls=":", label=r"LSA, any $\beta$")]
    h += [Line2D([], [], color=SCOL[k], marker="o", ls="", ms=3, label=f"DNS {SLAB[k]}") for k in ("1e-5", "1e-9", "0")]
    fig.legend(handles=h, loc="upper center", ncol=5, fontsize=6.5, bbox_to_anchor=(0.5, 1.07))
    fig.tight_layout(w_pad=0.4)
    save(fig, "sigma.png")
    for r in rows:
        print("w%d x-w=%d  LSAbox %.3f  LSAall %.3f  DNS1e-9 %.3f  DNS0 %.3f" % r)


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "dns":
    fig_dns()


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "sigma":
    fig_sigma()


# ---- 5. sigma(beta), fundamental and subharmonic ---------------------------------------------------
def fig_beta():
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.4), sharey=True)
    B, G = D["betas"], list(D["gammas"])
    XS = {50: "C0", 150: "C2", 300: "C3"}
    for j, w in enumerate((65, 66, 70)):
        iw = S.W.index(w)
        for X, c in XS.items():
            ist = list(ST).index(X)
            for g, ls in ((0.0, "-"), (0.5, "--")):
                s = D["lam"][iw, ist, G.index(g)].real / D["c"][iw, ist]
                ax[j].plot(B, s, color=c, ls=ls, marker="o" if g == 0 else None, ms=2)
                k = np.nanargmax(s)
                print(f"w{w} x-w={X} gamma={g}: beta_max={B[k]:.2f} sigma_max={s[k]:.3f}  box max={np.nanmax(s[S.BK]):.3f}")
        for k in range(1, 5):
            ax[j].axvline(2 * np.pi * k / 6, color="0.5", lw=0.6, ls=":")
        ax[j].axhline(0, color="0.4", lw=0.6)
        ax[j].set_title(rf"$w={w}$", fontsize=8.5); ax[j].set_xlabel(r"$\beta$"); ax[j].set_xlim(0, B.max())
    ax[0].set_ylabel(r"$\sigma=\Re\lambda/c$")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=c, label=rf"$x-w={X}$") for X, c in XS.items()]
    h += [Line2D([], [], color="k", marker="o", ms=2, label=r"$\gamma=0$ (fundamental)"),
          Line2D([], [], color="k", ls="--", label=r"$\gamma=1/2$ (subharmonic)"),
          Line2D([], [], color="0.5", ls=":", label=r"box $\beta_k=2\pi k/6$")]
    fig.legend(handles=h, loc="upper center", ncol=6, fontsize=6.3, bbox_to_anchor=(0.5, 1.08))
    fig.tight_layout(w_pad=0.4)
    save(fig, "beta.png")


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "beta":
    fig_beta()
