"""
N factor of a gap flow predicted from local (parallel flow) stability theory.

The model
--------
The perturbation at a station is written as an incoherent superposition of
Tollmien-Schlichting modes of *real* frequency,

    u'(x, y, t) = sum_j a_j exp(N_j(x)) phi_j(x, y) exp(-i omega_j t) + c.c.

with

    N_j(x) = - int_{x0}^{x} alpha_i(xi; omega_j) dxi

the usual e^N exponent of the mode of frequency omega_j, ``phi_j(x, .)`` the
local Orr-Sommerfeld eigenfunction at that station *normalised to unit L2 norm*
and ``a_j`` a weight that is constant in x.

Normalising the eigenfunction at every station is what makes the model
closed: the shape of the mode then carries no amplitude information, all of it
sits in exp(N_j), and everything the local theory cannot know -- receptivity,
the efficiency with which the forcing excites each frequency, the non-parallel
and non-modal corrections near the gap -- is pushed into the single set of
weights a_j.

Why the modes add in energy and not in amplitude
------------------------------------------------
The DNS amplitude this is compared against is

    A(x)^2 = int <u'^2 + v'^2> dy,

a *time averaged* quantity (``pp.DeltaN_computation.computeAmplitude``). Modes
of distinct real frequency are uncorrelated in time, so every cross term
averages to zero and, with ||phi_j|| = 1,

    A(x)^2 = sum_j a_j^2 exp(2 N_j(x)),
    N(x)   = 1/2 log( A(x)^2 / A(x0)^2 ).

So the eigenfunctions drop out of the prediction entirely and only alpha_i
enters -- which is the whole point of the unit norm convention, and is also why
the result is insensitive to the quadrature used for the norm.

Conventions
-----------
Everything is in the global non dimensionalisation of the DNS (lengths in
delta*_le, velocity in U_inf, ``re`` = Re_delta*_le), so alpha and omega need
no rescaling between stations and N_j is a plain integral in x. Frequencies are
labelled by the reduced frequency

    F = omega * 10^6 / re,

which is the invariant of a mode travelling downstream (omega and Re both scale
with the local delta*, so their ratio does not).
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from pp.colors import colors
from pp.orrSommerfeld import OrrSommerfeld

__all__ = [
    "Station",
    "load_sections",
    "blasius_stations",
    "track_modes",
    "ModeSet",
    "synthesise",
    "fit_weights",
]

BLASIUS_C = 1.7207876573


# --------------------------------------------------------------------------
# base flow stations
# --------------------------------------------------------------------------
@dataclass
class Station:
    x: float
    y: np.ndarray
    u: np.ndarray


def load_sections(
    data_dir: str,
    pattern: str = "points_n600_x*.dat",
    xmin: float | None = None,
    xmax: float | None = None,
) -> list[Station]:
    """
    Read the base flow profiles written by ``createPointsOfSectionDomain*``.

    The files are ``x y u v p`` with three header lines; ``y`` starts at the
    gap floor for a station over the gap and at the wall elsewhere.
    """
    files = sorted(Path(data_dir).glob(pattern))
    out: list[Station] = []
    for f in files:
        tag = f.stem.split("_x")[-1]
        try:
            x = float(tag)
        except ValueError:
            print(colors.WARNING + f"skipping unparsable {f.name}" + colors.ENDC)
            continue
        if xmin is not None and x < xmin:
            continue
        if xmax is not None and x > xmax:
            continue
        d = np.loadtxt(f, skiprows=3)
        out.append(Station(x, d[:, 1], d[:, 2]))
    out.sort(key=lambda s: s.x)
    if not out:
        raise FileNotFoundError(f"no sections matching {pattern} in {data_dir}")
    return out


def blasius_stations(
    x: np.ndarray,
    re: float,
    blasius_file: str = "/home/victor/Desktop/orrSommerfeldSolver/data/blasius.dat",
) -> list[Station]:
    """
    Analytical flat plate reference: the similarity profile stretched by the
    local displacement thickness delta*(x)/delta*_le = sqrt(1 + x C^2 / Re).

    Used for the flat plate baseline when no DNS sections were extracted for
    it, so that Delta N can still be formed against a like for like local
    theory curve.
    """
    d = np.loadtxt(blasius_file, skiprows=2)
    eta, u = d[:, 0], d[:, 1]  # profile already has delta* = 1
    out = []
    for xi in np.atleast_1d(x).astype(float):
        r = np.sqrt(1.0 + xi * BLASIUS_C**2 / re)
        out.append(Station(float(xi), eta * r, u.copy()))
    return out


# --------------------------------------------------------------------------
# mode tracking
# --------------------------------------------------------------------------
@dataclass
class ModeSet:
    """alpha(x, omega) for a set of frequencies, plus the resulting N_j(x)."""

    x: np.ndarray  # (nx,)
    F: np.ndarray  # (nf,) reduced frequency
    omega: np.ndarray  # (nf,) global omega
    alpha: np.ndarray  # (nf, nx) complex, nan where the mode was lost
    n_global: int = 0  # number of (station, frequency) eigensolves performed
    meta: dict = field(default_factory=dict)

    @property
    def growth(self) -> np.ndarray:
        """-alpha_i, positive meaning amplified."""
        return -self.alpha.imag

    def nfactor(
        self,
        clip_stable_upstream: bool = False,
        freeze: tuple[float, float] | None = None,
    ) -> np.ndarray:
        """
        N_j(x) = int -alpha_i dx, cumulative from the first station.

        ``clip_stable_upstream`` starts each mode's integral at its own branch
        I instead, i.e. the textbook envelope convention. It is *off* by
        default: the DNS reference amplitude is taken at a fixed upstream
        station and includes the decay of the damped modes, so integrating from
        the same place for every mode is the consistent choice.

        ``freeze=(x0, x1)`` drops the contribution of that interval, leaving
        N_j flat across it. Passing the gap is the point: a parallel flow
        analysis of a station *inside* a short gap is the weakest link in the
        whole model -- the shear layer there is a few gap widths long, far from
        the slowly varying limit the local theory assumes -- so being able to
        run with and without it separates "the gap amplifies" from "the local
        theory says the gap amplifies by this much".
        """
        g = np.nan_to_num(self.growth, nan=0.0)
        if freeze is not None:
            inside = (self.x >= freeze[0]) & (self.x <= freeze[1])
            g = np.where(inside[None, :], 0.0, g)
        n = np.zeros_like(g)
        n[:, 1:] = np.cumsum(0.5 * (g[:, 1:] + g[:, :-1]) * np.diff(self.x), axis=1)
        if clip_stable_upstream:
            for j in range(n.shape[0]):
                pos = np.flatnonzero(g[j] > 0.0)
                n[j] -= n[j, pos[0]] if pos.size else n[j, -1]
                n[j] = np.maximum(n[j], 0.0)
        return n


# selection window for the discrete instability mode. The phase speed is what
# separates it from everything else: a Tollmien-Schlichting mode runs at
# c_r = 0.2..0.5 and the shear layer mode over the gap at c_r = 0.3..0.5, while
# the discretised continuous spectrum piles up at c_r -> 1.
C_MIN, C_MAX = 0.15, 0.75
# |alpha_i| < RATIO * alpha_r, i.e. no more than ~1 e-fold per RATIO*2pi
# wavelengths. Rejects the strongly damped junk the companion pencil produces
# without cutting into the real gap modes, which reach |alpha_i|/alpha_r ~ 0.6.
RATIO_MAX = 0.9
# Backstop on the wavenumber, as a multiple of omega / C_MIN. It has to be
# written relative to omega and not as an absolute number: omega = F * re * 1e-6
# grows with the Reynolds number, so alpha_r ~ omega / c_r does too, and a fixed
# cap that is loose at re = 1000 (alpha_r reaches 2) silently deletes the real
# mode at re = 3000 (alpha_r reaches 4.3 at the top of the band). Worse than
# deleting it: select_mode would then return the least damped *other* admissible
# eigenvalue rather than nothing, so the error would not announce itself.
ALPHA_MAX_FACTOR = 1.5


def select_mode(alpha: np.ndarray, omega: float) -> int | None:
    """
    Index of the physical instability mode in a full spatial spectrum.

    The spatial spectrum holds, besides the mode we want, the upstream running
    branch, the discretised continuous spectrum and the spurious eigenvalues
    any companion linearisation produces. Filtering on the phase speed and on
    the damping per wavelength leaves the discrete branch, and the least damped
    survivor is the one that sets the N factor.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        c = omega / alpha.real
    ok = (
        (alpha.real > 1e-6)
        & (np.abs(alpha) < ALPHA_MAX_FACTOR * omega / C_MIN)
        & np.isfinite(c)
        & (c > C_MIN)
        & (c < C_MAX)
        & (np.abs(alpha.imag) < RATIO_MAX * alpha.real)
    )
    idx = np.flatnonzero(ok)
    if idx.size == 0:
        return None
    return int(idx[np.argmax(-alpha.imag[idx])])


def _scan_station(
    args: tuple[Station, float, np.ndarray, int, float],
) -> tuple[np.ndarray, float]:
    """One station, every frequency. Runs in a worker process."""
    st, re, omega, n, ymax = args
    os_ = OrrSommerfeld(st.y, st.u, re=re, n=n, ymax=ymax)
    out = np.full(omega.size, np.nan + 1j * np.nan, dtype=complex)
    for j, om in enumerate(omega):
        al, _ = os_.spatial(om, vectors=False)
        k = select_mode(al, om)
        if k is not None:
            out[j] = al[k]
    return out, os_.fit_rms


def track_modes(
    stations: list[Station],
    re: float,
    F: np.ndarray,
    n: int = 101,
    ymax: float = 45.0,
    workers: int | None = None,
    verbose: bool = True,
) -> ModeSet:
    """
    Solve the spatial problem at every station and every frequency.

    Each (station, frequency) is solved *independently* with a full
    eigendecomposition and the mode is chosen by :func:`select_mode`. Marching
    a single eigenvalue with a Newton continuation would be about a hundred
    times cheaper, but it assumes the branch it started on stays the relevant
    one -- which is exactly what fails over the gap, where the shear layer mode
    appears and overtakes the Tollmien-Schlichting one. Solving independently
    costs a few minutes per case and makes no such assumption; the stations are
    embarrassingly parallel.
    """
    F = np.atleast_1d(np.asarray(F, float))
    omega = F * re * 1e-6
    x = np.array([s.x for s in stations], float)

    args = [(st, float(re), omega, n, ymax) for st in stations]
    if workers is None:
        workers = min(len(stations), os.cpu_count() or 1)

    if workers > 1:
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(max_workers=workers) as ex:
            res = list(ex.map(_scan_station, args, chunksize=1))
    else:
        res = [_scan_station(a) for a in args]

    alpha = np.array([r[0] for r in res]).T  # (nf, nx)
    fit_rms = np.array([r[1] for r in res])

    lost = int(np.sum(~np.isfinite(alpha)))
    if verbose and lost:
        print(
            colors.WARNING
            + f"  no admissible mode at {lost}/{alpha.size} (station, frequency) "
            "pairs; those are treated as neutral in the integral"
            + colors.ENDC
        )

    return ModeSet(
        x=x,
        F=F,
        omega=omega,
        alpha=alpha,
        n_global=alpha.size,
        meta={
            "re": re,
            "n": n,
            "ymax": ymax,
            "fit_rms_max": float(np.max(fit_rms)),
            "lost": lost,
        },
    )


def jump_diagnostic(ms: ModeSet, rel: float = 0.35) -> list[tuple[float, float]]:
    """
    (F, x) pairs where alpha changes by more than ``rel`` of its own size
    between neighbouring stations.

    Because every station is solved independently, a jump is a genuine signal
    that the dominant branch has changed -- worth seeing rather than smoothing
    away, since that is what the gap does to the spectrum.
    """
    out = []
    a = ms.alpha
    for j in range(a.shape[0]):
        d = np.abs(np.diff(a[j]))
        scale = np.maximum(np.abs(a[j, :-1]), 1e-9)
        for i in np.flatnonzero(d > rel * scale):
            out.append((float(ms.F[j]), float(ms.x[i + 1])))
    return out


# --------------------------------------------------------------------------
# synthesis
# --------------------------------------------------------------------------
def synthesise(
    ms: ModeSet,
    weights: np.ndarray | None = None,
    clip_stable_upstream: bool = False,
    freeze: tuple[float, float] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Combine the per frequency N_j into the N factor of the summed field.

        A(x)^2 = sum_j a_j^2 exp(2 N_j(x)),
        N(x)   = 1/2 log( A^2(x) / A^2(x0) ),

    with x0 chosen exactly as ``pp.DeltaN_computation.computeNx`` chooses it:
    the station of smallest amplitude upstream of the gap. Returns ``(x, N)``
    truncated to x >= x0 so that the two curves can be compared directly.

    ``weights=None`` gives every frequency the same weight, which is the
    zeroth order model for the white in time forcing the DNS uses.
    """
    N = ms.nfactor(clip_stable_upstream=clip_stable_upstream, freeze=freeze)
    a = np.ones(ms.F.size) if weights is None else np.asarray(weights, float)
    if a.shape != ms.F.shape:
        raise ValueError(f"weights must have shape {ms.F.shape}, got {a.shape}")

    # log-sum-exp, because N_j reaches O(10) and exp(2N) overflows otherwise
    e = 2.0 * N + 2.0 * np.log(np.maximum(a, 1e-300))[:, None]
    m = np.max(e, axis=0)
    logA2 = m + np.log(np.sum(np.exp(e - m), axis=0))

    upstream = np.flatnonzero(ms.x <= 0.0)
    i0 = int(upstream[np.argmin(logA2[upstream])]) if upstream.size else 0
    return ms.x[i0:], 0.5 * (logA2[i0:] - logA2[i0])


def fit_weights(
    ms: ModeSet,
    x_dns: np.ndarray,
    n_dns: np.ndarray,
    clip_stable_upstream: bool = False,
    xfit: tuple[float, float] | None = None,
) -> np.ndarray:
    """
    Least squares weights a_j^2 >= 0 that make the synthesised amplitude follow
    the DNS one.

    A(x)^2 is *linear* in the unknowns a_j^2, so fitting A^2 rather than N is
    a non negative least squares problem with no local minima -- the price is
    that it weights the downstream, large amplitude end of the curve most,
    which is the part Delta N is read off anyway.
    """
    from scipy.optimize import nnls

    N = ms.nfactor(clip_stable_upstream=clip_stable_upstream)
    x = ms.x
    m = np.ones(x.size, bool) if xfit is None else (x >= xfit[0]) & (x <= xfit[1])

    n_i = np.interp(x, x_dns, n_dns)
    # work in units of the DNS amplitude so the fit is not dominated by scale
    lg = 2.0 * N[:, m]
    shift = np.max(lg, axis=0)
    design = np.exp(lg - shift)  # (nf, nm)
    target = np.exp(2.0 * n_i[m] - shift)

    q, _ = nnls(design.T, target)
    return np.sqrt(q)
