"""Herbert (1988)-style plots of the amplitude sweep of the secondary instability.

Diagnostic for docs/notes/20261108_secondaryStability.md, from
data/stability/secondaryStabilityAmplitude_d1.5.npz (scripts/secondaryStabilityAmplitude.py).
Run from scripts/pp: .venv/bin/python plots/secondaryAmplitude.py [X_SIGMA=100]
Writes to docs/notes/img/secondaryStability/:
  herbert_sigma.png      sigma(beta) for several wave amplitudes a at x - w = X_SIGMA
  herbert_neutral.png    neutral curves sigma = 0 in the (beta, a) plane at each station
  herbert_amplitude.png  sigma_max(a) for all widths and stations, and a_DNS(x)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = np.load(os.path.join(HERE, "../../../data/stability/secondaryStabilityAmplitude_d1.5.npz"))
NF = np.load(os.path.join(HERE, "../../../data/stability/secondaryStabilityNfactor_d1.5.npz"))
OUT = os.path.join(HERE, "../../../docs/notes/img/secondaryStability")
XS = int(sys.argv[1]) if len(sys.argv) > 1 else 100
W, ST, B, G, A = D["widths"], D["stations"], D["betas"], D["gammas"], D["amps"]
A0 = 1e-3                       # where a = 0 is drawn on log axes
Aplot = np.where(A > 0, A, A0)
LS = {0.0: "-", 0.5: "--"}
COL = {60: "#9ecae1", 63: "#4292c6", 64: "#08519c", 65: "#fdae6b", 66: "#e6550d", 70: "#7f2704"}
# stations at and downstream of these x - w are dropped: the DNS wave has died and the fit is noise
DEAD = {60: 200, 63: 300, 64: 450}
os.makedirs(OUT, exist_ok=True)


def valid(iw, ix):
    return ST[ix] < DEAD.get(int(W[iw]), 1e9) and not np.all(np.isnan(D["sigma"][iw, ix]))


def neutral(s):
    """Neutral amplitude along a (s: sigma over AMPS, nan = skipped): 0 if unstable at a = 0,
    (largest a tested, True) if stable at every a tested."""
    m = np.isfinite(s)
    s, amps = s[m], A[m]
    if not len(s):
        return np.nan, False
    if s[0] > 0:
        return 0.0, False
    for i in range(1, len(s)):
        if s[i] > 0:
            a0, a1 = max(amps[i - 1], A0), amps[i]
            t = -s[i - 1] / (s[i] - s[i - 1])
            return float(np.exp(np.log(a0) + t * np.log(a1 / a0))), False
    return float(amps[-1]), True


# ---- 1. sigma(beta) at x - w = XS for several a ----------------------------------------
ix = int(np.argmin(np.abs(ST - XS)))
cm = plt.get_cmap("viridis")
fig, axs = plt.subplots(2, 3, figsize=(15, 7.5), sharex=True)
for iw, (w, ax) in enumerate(zip(W, axs.flat)):
    if not valid(iw, ix):
        continue
    for ia, a in enumerate(A):
        for ig, g in enumerate(G):
            s = D["sigma"][iw, ix, ig, ia]
            if np.all(np.isnan(s)):
                continue
            ax.plot(B, s, LS[g], color=cm(ia / (len(A) - 1)), lw=1, label=f"a={a:g}" if ig == 0 else None)
    for ig, g in enumerate(G):
        ax.plot(B, D["sigma_dns"][iw, ix, ig], LS[g], color="k", lw=2.2,
                label=f"DNS a={D['a_dns'][iw, ix]:.3f}" if ig == 0 else None)
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_title(f"w={w}, x-w={ST[ix]}  (solid fundamental, dashed subharmonic)", fontsize=9)
    ax.legend(fontsize=6.5, ncol=2)
    ax.set_ylim(-0.04, 0.13)
for ax in axs[-1]:
    ax.set_xlabel("spanwise wavenumber beta")
for ax in axs[:, 0]:
    ax.set_ylabel("spatial growth rate sigma = Re(sigma_t)/c")
fig.tight_layout(); fig.savefig(f"{OUT}/herbert_sigma.png", dpi=110); plt.close(fig)

# ---- 2. neutral curves in (beta, a) -------------------------------------------------------
fig, axs = plt.subplots(2, 3, figsize=(15, 7.5), sharex=True, sharey=True)
cs = plt.get_cmap("plasma")
for iw, (w, ax) in enumerate(zip(W, axs.flat)):
    for ix, x0 in enumerate(ST):
        if not valid(iw, ix):
            continue
        col = cs(ix / len(ST))
        for ig, g in enumerate(G):
            r = [neutral(D["sigma"][iw, ix, ig, :, ib]) for ib in range(len(B))]
            a = np.array([max(v[0], A0) for v in r]); lb = np.array([v[1] for v in r])
            ax.plot(B, a, LS[g], color=col, lw=1.2, label=f"x-w={x0}" if ig == 0 else None)
            ax.plot(B[lb], a[lb], "^", color=col, ms=4)
        ax.axhline(D["a_dns"][iw, ix], color=col, lw=0.6, ls=":")
    ax.set_yscale("log"); ax.set_xscale("log")
    ax.set_title(f"w={w}: neutral curves sigma=0, unstable above (solid fund., dashed subh.)\n"
                 f"dotted: DNS amplitude at that station", fontsize=8)
    ax.legend(fontsize=6.5)
ax.set_ylim(A0 * 0.9, 0.5)
for ax in axs[-1]:
    ax.set_xlabel("beta")
for ax in axs[:, 0]:
    ax.set_ylabel("wave amplitude a = max 2|u_1|   (a = 0 drawn at 1e-3)")
fig.tight_layout(); fig.savefig(f"{OUT}/herbert_neutral.png", dpi=110); plt.close(fig)

# ---- 3. sigma_max(a) for every width and station, and a_DNS(x) ------------------------
MK = dict(zip(ST, "osD^vP"))
fig, axs = plt.subplots(1, 3, figsize=(19, 5))
ib1 = int(np.argmin(np.abs(B - 2 * np.pi / 6)))
for iw, w in enumerate(W):
    c = COL[int(w)]
    for ix, x0 in enumerate(ST):
        if not valid(iw, ix):
            continue
        smax = np.nanmax(D["sigma"][iw, ix], axis=(0, 2))
        s1 = np.nanmax(D["sigma"][iw, ix, :, :, ib1], axis=0)
        m = np.isfinite(smax)
        for ax, s, sd in ((axs[0], smax, np.nanmax(D["sigma_dns"][iw, ix])),
                          (axs[1], s1, np.nanmax(D["sigma_dns"][iw, ix, :, ib1]))):
            ax.plot(Aplot[m], s[m], "-", marker=MK[x0], color=c, ms=4, lw=0.7)
            ax.plot(D["a_dns"][iw, ix], sd, marker=MK[x0], color=c, ms=10, mec="k", ls="none")
    jw = list(NF["widths"]).index(w)
    am = NF["amp"][jw].copy()
    if int(w) in DEAD:
        am[NF["stations"] >= DEAD[int(w)] - 50] = np.nan
    axs[2].plot(NF["stations"], am, color=c, lw=1.8, label=f"w={w}")
for ax, t in zip(axs[:2], ("max over beta and gamma (fundamental and subharmonic)",
                           f"beta = {B[ib1]:.2f} (grid value nearest 2 pi/6 = 1.05), max over gamma")):
    ax.set_xscale("log"); ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xlabel("2D wave amplitude a = max 2|u_1|  (a = 0 drawn at 1e-3)")
    ax.set_title(f"secondary growth rate against amplitude, {t}\n"
                 "small: rescaled wave, one line per (w, station); large: the DNS wave", fontsize=8.5)
    ax.set_xlim(8e-4, 0.4)
axs[0].set_ylabel("spatial growth rate sigma")
for x0, mk in MK.items():
    axs[0].plot([], [], mk, color="0.5", label=f"x-w={x0}")
for w in W:
    axs[0].plot([], [], "-", color=COL[int(w)], label=f"w={w}")
axs[0].legend(fontsize=7, ncol=2)
for lev in (0.02, 0.04):    # amplitude giving sigma_max = lev, median over widths and stations
    av = []
    for iw in range(len(W)):
        for ix in range(len(ST)):
            if not valid(iw, ix):
                continue
            smax = np.nanmax(D["sigma"][iw, ix], axis=(0, 2))
            m = np.isfinite(smax) & (A > 0)
            if smax[m].max() > lev:
                av.append(np.exp(np.interp(lev, smax[m], np.log(A[m]))))
    am = np.median(av)
    print(f"sigma_max = {lev}: a = {am:.4f} (range {min(av):.4f}-{max(av):.4f})")
    axs[2].axhline(am, color="0.4", ls=":", lw=0.8)
    axs[2].text(5, am, f"sigma_max = {lev}", va="bottom", fontsize=7)
axs[2].set_yscale("log"); axs[2].set_ylim(1e-3, 0.4); axs[2].set_xlim(0, 600)
axs[2].set_xlabel("x - w"); axs[2].set_ylabel("a_DNS(x) = max 2|u_1| of the 2D limit cycle")
axs[2].set_title("2D wave amplitude along x (cut where the wave has died)", fontsize=9)
axs[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{OUT}/herbert_amplitude.png", dpi=110); plt.close(fig)
print("saved", OUT)
