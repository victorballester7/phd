from pathlib import Path
import os
import subprocess
import warnings

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.axes import Axes

from pp.fileManagement import editFile, extract_depth_width

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

BLASIUS_C = 1.7207876573

ALPHA_MIN = 0.01
ALPHA_MAX = 0.5
ALPHA_NUM = 45

OMEGA_MIN = 0.01
OMEGA_MAX = 0.5
OMEGA_NUM = 30

# Order of the Gaster transformation (1 or 2). Order 2 keeps the curvature of
# the temporal dispersion relation and stays usable when the amplification is
# no longer small; order 1 is the classical alpha_i = -omega_i / c_g.
GASTER_ORDER = 2

# Smoothing applied to omega_r(alpha_r), omega_i(alpha_r) before differentiating.
# None differentiates the raw sweep; a float is the `lam` of a smoothing spline
# and is worth setting when the eigenvalue sweep is ragged, since the order 2
# transformation needs second derivatives.
GASTER_SPLINE_LAM: float | None = None

# What to do when the output file of a run already exists:
#   "resume"    reuse every station already stored and only run the solver for
#               the ones that are missing (default)
#   "reuse"     plot what is stored and never call the solver
#   "overwrite" recompute everything
#   "ask"       prompt, as the script used to
CACHE_POLICY = "resume"

# absolute tolerance when matching a station against the cached ones
X_TOL = 1e-6


def run_os(toml_file: str) -> None:
    solver_dir = Path(toml_file).parent.parent
    subprocess.run(["make", "run"], cwd=solver_dir, check=True)


def setup_toml(
    filename_toml: str, alpha_r_min: float, alpha_r_max: float, alpha_r_num: int
) -> None:
    line_startswith = np.array(
        [
            "n = ",
            "re = ",
            "beta = ",
            "useTargetEV = ",
            "vars_r = ",
            "vars_i = ",
            "branch = ",
            "doPlot = ",
            "use_c = ",
            "multipleRun = ",
            "plotUprofile = ",
            "colX = ",
            "colY = ",
            "numSkipHeaderLines = ",
        ]
    )
    replacement_line = np.array(
        [
            "n = 175",
            "re = 1000",
            "beta = { r = 0.0, i = 0.0 }",
            "useTargetEV = false",
            f"vars_r = {{min = {alpha_r_min}, max = {alpha_r_max}, num = {alpha_r_num}}}",
            "vars_i = {min = 0.0, max = 0.0, num = 1}",
            'branch = "spatial"',
            "doPlot = false",
            "use_c = false",
            "multipleRun = true",
            "plotUprofile = false",
            "colX = 1",
            "colY = 2",
            "numSkipHeaderLines = 1",
        ]
    )
    editFile(filename_toml, replacement_line, line_startswith)


def read_evs(
    filename: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read temporal OS output assuming columns: ev_r, ev_i, var_r, ...
    """
    data = np.loadtxt(filename, skiprows=1)
    data = np.atleast_2d(data)
    if data.shape[1] < 3:
        raise ValueError(
            f"Expected at least 3 columns in {filename}, got {data.shape[1]}"
        )
    ev_r = data[:, 0]
    ev_i = data[:, 1]
    var_r = data[:, 2]
    return ev_r, ev_i, var_r


def _interpolate_zero_crossing(v1: float, y1: float, v2: float, y2: float) -> float:
    return (v1 * y2 - v2 * y1) / (y2 - y1)


def get_neutral_var_r(
    ev_i: np.ndarray,
    var_r_solver: np.ndarray,
) -> np.ndarray:
    """
    Convert neutral condition from c_i = 0 to omega_r.

    In temporal analysis with real alpha_r:
        omega = alpha * c
        omega_i = alpha_r * c_i
    so omega_i = 0 is equivalent to c_i = 0 (for alpha_r > 0).
    """
    var_r_solver_neutralCurve: list[float] = []
    nvals = len(ev_i)
    if nvals < 2:
        return np.array([], dtype=float)

    for idx in range(nvals - 1):
        ev_i1 = ev_i[idx]
        ev_i2 = ev_i[idx + 1]

        if ev_i1 == 0.0:
            var_r_solver_neutralCurve.append(var_r_solver[idx])
            continue

        if ev_i1 * ev_i2 < 0.0:
            ev_r_neutral = _interpolate_zero_crossing(
                var_r_solver[idx], ev_i1, var_r_solver[idx + 1], ev_i2
            )
            var_r_solver_neutralCurve.append(ev_r_neutral)

    return np.asarray(var_r_solver_neutralCurve, dtype=float)




def _second_derivative(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """
    Second derivative on a possibly non-uniform grid, exact for quadratics.

    Differentiating ``np.gradient`` twice would instead give a 2h-wide stencil
    that both loses accuracy and leaves noise untouched, which matters here
    because the order 2 Gaster transformation leans on this quantity.
    """
    h1 = x[1:-1] - x[:-2]
    h2 = x[2:] - x[1:-1]

    interior = (
        2.0
        * (h2 * y[:-2] - (h1 + h2) * y[1:-1] + h1 * y[2:])
        / (h1 * h2 * (h1 + h2))
    )

    # one sided formulas are too noisy at the ends of the sweep to be worth it
    return np.concatenate(([interior[0]], interior, [interior[-1]]))


def _smoothed(y: np.ndarray, x: np.ndarray, lam: float | None) -> np.ndarray:
    if lam is None or len(x) < 5:
        return y
    from scipy.interpolate import make_smoothing_spline

    return make_smoothing_spline(x, y, lam=lam)(x)


def gasterTransform(
    omega_r: np.ndarray,
    omega_i: np.ndarray,
    alpha_r: np.ndarray,
    order: int = GASTER_ORDER,
    lam: float | None = GASTER_SPLINE_LAM,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Map a temporal sweep onto the spatial branch.

    The temporal solver gives the analytic function ω(α) sampled on the real
    axis, α = a. The spatial mode at that station is α = a + i b with ω real,
    so expanding about the real axis,

        ω(a + ib) = ω(a) + i b ω'(a) - (b²/2) ω''(a) + O(b³)

    and splitting into real and imaginary parts,

        Im:  ω_i + b ω_r' - (b²/2) ω_i'' = 0          -> b = α_i
        Re:  Ω  = ω_r - b ω_i' - (b²/2) ω_r''         -> real frequency

    *order = 1* drops every quadratic term and the ``b ω_i'`` frequency shift,
    recovering the classical α_i = -ω_i / c_g at fixed ω = ω_r. It is only
    justified while the amplification is small, since the neglected terms grow
    like b².

    *order = 2* solves the imaginary equation as the quadratic it is,

        b = -2 ω_i / (ω_r' + sqrt(ω_r'² + 2 ω_i'' ω_i))

    written in the form that avoids the cancellation of the textbook root when
    ω_i'' is small, and reduces to the order 1 result exactly at ω_i'' = 0.
    The frequency is corrected to the matching order. Where the discriminant
    turns negative the expansion has no real spatial branch and that point
    falls back to order 1.

    Returns ``(alpha_i, omega, alpha_i_first_order)``; comparing the last two
    shows how hard the transformation is working.
    """
    if order not in (1, 2):
        raise ValueError(f"Gaster order must be 1 or 2, got {order}")

    omega_r_s = _smoothed(omega_r, alpha_r, lam)
    omega_i_s = _smoothed(omega_i, alpha_r, lam)

    # edge_order=2 matters: the default first order one sided difference leaves
    # an O(h) error in c_g at the two ends of the swept range, which is where
    # the most amplified mode often sits
    cg = np.gradient(omega_r_s, alpha_r, edge_order=2)
    # |c_g| guards against the sign flips an occasional misconverged eigenvalue
    # introduces in the sweep.
    cg = np.abs(cg)

    with np.errstate(divide="ignore", invalid="ignore"):
        alpha_i_first = -omega_i_s / cg

    if order == 1:
        return alpha_i_first, omega_r.copy(), alpha_i_first

    domega_i = np.gradient(omega_i_s, alpha_r, edge_order=2)
    d2omega_i = _second_derivative(omega_i_s, alpha_r)
    d2omega_r = _second_derivative(omega_r_s, alpha_r)

    disc = cg**2 + 2.0 * d2omega_i * omega_i_s
    denom = cg + np.sqrt(np.where(disc > 0.0, disc, 0.0))

    usable = (disc > 0.0) & (denom != 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        alpha_i = np.where(usable, -2.0 * omega_i_s / denom, alpha_i_first)

    n_fallback = int(np.sum(~usable & np.isfinite(alpha_i_first)))
    if n_fallback:
        print(
            f"  Gaster order 2: {n_fallback}/{len(alpha_r)} points had no real "
            "spatial branch, order 1 used there."
        )

    omega = omega_r_s - alpha_i * domega_i - 0.5 * alpha_i**2 * d2omega_r
    omega = np.where(usable, omega, omega_r)

    return alpha_i, omega, alpha_i_first

def _refine_max(s: np.ndarray, g: np.ndarray) -> tuple[float, float]:
    """
    Parabolic refinement of the maximum of the sampled curve g(s).

    Falls back to the sampled point itself when the discrete maximum sits on
    the edge of the swept range (the true maximum is then outside it), or when
    the three-point fit is not concave.
    """
    i = int(np.argmax(g))

    if i == 0 or i == len(g) - 1:
        return float(s[i]), float(g[i])

    coeffs = np.polyfit(s[i - 1 : i + 2], g[i - 1 : i + 2], 2)

    if coeffs[0] >= 0.0:
        return float(s[i]), float(g[i])

    s_star = -coeffs[1] / (2.0 * coeffs[0])
    if not (s[i - 1] <= s_star <= s[i + 1]):
        return float(s[i]), float(g[i])

    return float(s_star), float(np.polyval(coeffs, s_star))


def get_maxGrowthRate(
    ev_r: np.ndarray,
    ev_i: np.ndarray,
    var_r: np.ndarray,
    temporal: bool,
    useGaster: bool,
    gasterOrder: int = GASTER_ORDER,
    gasterLam: float | None = GASTER_SPLINE_LAM,
) -> tuple[float, float, float] | None:
    """
    Locate the most amplified mode of a single station.

    Returns ``(var_r*, ev_r*, growth*)`` with ``growth*`` always the physical
    growth rate, positive meaning amplified:

    * temporal            : sweep α_r, growth = ω_i    -> (α_r*, ω_r*, ω_i*)
    * temporal + Gaster   : sweep α_r, growth = -α_i = ω_i / c_g, reported as a
                            spatial result parametrised by ω_r
                                               -> (ω_r*, α_r*, -α_i*)
    * spatial             : sweep ω_r, growth = -α_i   -> (ω_r*, α_r*, -α_i*)

    The sign convention matches ``plot_max_growthRate`` / ``saveData``, whose
    first column is the swept variable of the *reported* branch.

    ``None`` is returned when the station has too few usable eigenvalues. A
    station whose maximum is still damped is *not* discarded: it is returned
    with a negative ``growth*`` so that the locus stays continuous through the
    stable part of the domain.
    """
    order = np.argsort(var_r)
    var_r, ev_r, ev_i = var_r[order], ev_r[order], ev_i[order]

    if temporal and useGaster:
        alpha_i, omega, alpha_i_first = gasterTransform(
            ev_r, ev_i, var_r, gasterOrder, gasterLam
        )
        growth = -alpha_i

        good = np.isfinite(growth) & np.isfinite(omega)
        if good.sum() < 3:
            return None

        alpha_r = var_r[good]
        alpha_r_star, growth_star = _refine_max(alpha_r, growth[good])
        omega_star = float(np.interp(alpha_r_star, alpha_r, omega[good]))

        if gasterOrder == 2:
            # how far the order 2 result has moved from order 1 at the peak:
            # a correction that is not small means the expansion itself is
            # marginal and the spatial branch should be solved directly
            first = float(np.interp(alpha_r_star, alpha_r, -alpha_i_first[good]))
            if abs(growth_star) > 0.0:
                rel = abs(growth_star - first) / abs(growth_star)
                if rel > 0.2:
                    print(
                        f"  Gaster order 2 correction is {100 * rel:.0f}% of the "
                        "growth rate; the expansion is marginal here."
                    )

        return omega_star, alpha_r_star, growth_star

    # Plain temporal (growth = omega_i) or plain spatial (growth = -alpha_i).
    growth = ev_i if temporal else -ev_i

    good = np.isfinite(growth)
    if good.sum() < 3:
        return None

    var_r, ev_r, growth = var_r[good], ev_r[good], growth[good]
    var_r_star, growth_star = _refine_max(var_r, growth)
    ev_r_star = float(np.interp(var_r_star, var_r, ev_r))

    return var_r_star, ev_r_star, growth_star



def _load_neutral_curve(path: str | None) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Load a neutral curve file and split it into its two branches.

    The files hold columns (Re_δ*/1.72, ω_r, α_r, x) with both neutral points
    of a station listed one after the other, so simply sorting by Re_δ* would
    zig-zag between branch I and branch II. Stations are grouped by x and the
    lower/upper neutral point of each is sent to its own branch, ordered by x.

    Returns ``(lower, upper)``, each with the same four columns.
    """
    if not path or not os.path.isfile(path):
        return None

    data = np.atleast_2d(np.loadtxt(path, skiprows=1))
    if data.shape[1] < 3:
        return None
    if data.shape[1] == 3:
        # older files were written without the alpha_r column
        data = np.column_stack(
            (data[:, 0], data[:, 1], np.full(len(data), np.nan), data[:, 2])
        )

    lower: list[np.ndarray] = []
    upper: list[np.ndarray] = []
    for x in np.unique(data[:, 3]):
        station = data[data[:, 3] == x]
        # branch membership is decided on omega_r and carried over to alpha_r,
        # so the two columns stay on the same branch
        order = np.argsort(station[:, 1])
        lower.append(station[order[0]])
        upper.append(station[order[-1]])

    return np.asarray(lower), np.asarray(upper)


def _coloured_line(
    ax: Axes,
    x: np.ndarray,
    y: np.ndarray,
    c: np.ndarray,
    vmax: float,
    cmap: str = "coolwarm",
    linewidth: float = 2.0,
) -> LineCollection:
    """Draw x-y as a single line whose colour follows c."""
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, cmap=cmap, norm=plt.Normalize(-vmax, vmax))
    lc.set_array(0.5 * (c[:-1] + c[1:]))  # colour each segment by its midpoint
    lc.set_linewidth(linewidth)
    lc.set_zorder(2)
    ax.add_collection(lc)

    return lc


def plot_max_growthRate(
    reynolds: np.ndarray,
    var_r: np.ndarray,
    ev_r: np.ndarray,
    growthRate: np.ndarray,
    xpositions: np.ndarray,
    temporal: bool,
    neutral_curve_file: str | None = None,
) -> None:
    """
    Plot the locus of the most amplified mode, coloured by its growth rate.

    Two panels show the locus in the (Re_δ*, ω_r) and (Re_δ*, α_r) planes -- the
    line that runs inside the neutral curve -- and a third shows the growth rate
    against the streamwise station x.
    """
    if reynolds.size == 0 or var_r.size == 0:
        raise ValueError(
            "No maximum growth rate points were found. Check the swept range "
            "and the solver output."
        )

    # the locus is a curve along the domain, so it must be ordered by station
    order = np.argsort(xpositions)
    reynolds = reynolds[order]
    var_r = var_r[order]
    ev_r = ev_r[order]
    growthRate = growthRate[order]
    xpositions = xpositions[order]

    var_label = "α" if temporal else "ω"
    ev_label = "ω" if temporal else "α"
    growth_label = "ω_i" if temporal else "-α_i"

    reynolds = BLASIUS_C * reynolds  # rescale to Re_δ*

    vmax = float(np.max(np.abs(growthRate)))
    if vmax == 0.0:
        vmax = 1.0

    neutral = _load_neutral_curve(neutral_curve_file)

    # The neutral curve files store (ω_r, α_r), i.e. a spatial parametrisation,
    # so they only line up with the panels when the result is reported
    # spatially (plain spatial branch, or temporal + Gaster).
    for column, (y, label) in enumerate(
        ((var_r, var_label), (ev_r, ev_label)), start=1
    ):
        fig, a = plt.subplots()

        if neutral is not None and not temporal:
            for branch, lab in zip(neutral, ("neutral curve", None)):
                yn = branch[:, column]
                if np.isfinite(yn).any():
                    a.plot(
                        BLASIUS_C * branch[:, 0],
                        yn,
                        "k-",
                        linewidth=1,
                        zorder=1,
                        label=lab,
                    )
            a.legend()

        lc = _coloured_line(a, reynolds, y, growthRate, vmax)
        a.scatter(
            reynolds,
            y,
            c=growthRate,
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
            s=12,
            edgecolors="k",
            linewidths=0.3,
            zorder=3,
        )
        fig.colorbar(lc, ax=a, label=f"growth rate {growth_label}")

        a.set_xlabel("Re_δ*")
        a.set_ylabel(label + "_r")
        a.set_title(
            "Maximum growth rate locus " + ("(Temporal)" if temporal else "(Spatial)")
        )
        a.grid(zorder=0)
        a.autoscale()
        fig.tight_layout()

    fig, a = plt.subplots()
    a.axhline(0.0, color="k", linewidth=0.8, zorder=1)
    _coloured_line(a, xpositions, growthRate, growthRate, vmax)
    a.scatter(
        xpositions,
        growthRate,
        c=growthRate,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
        s=12,
        edgecolors="k",
        linewidths=0.3,
        zorder=3,
    )
    a.set_xlabel("x")
    a.set_ylabel(f"max growth rate {growth_label}")
    a.grid(zorder=0)
    a.autoscale()
    fig.tight_layout()

    plt.show()



def _format_meta(meta: dict) -> str:
    return "meta: " + " ".join(f"{k}={v}" for k, v in sorted(meta.items()))


def _parse_meta(filename: str) -> dict[str, str]:
    """Read back the `# meta:` line written next to the column header."""
    with open(filename) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            body = line.lstrip("#").strip()
            if body.startswith("meta:"):
                return dict(
                    tok.split("=", 1)
                    for tok in body[len("meta:") :].split()
                    if "=" in tok
                )
    return {}


def resolveCache(
    filename_save: str,
    x_wanted: np.ndarray,
    meta: dict,
    cachePolicy: str = CACHE_POLICY,
) -> tuple[np.ndarray | None, np.ndarray, bool]:
    """
    Decide which stations still have to be sent to the solver.

    Every station costs one full Orr-Sommerfeld run, so an existing output file
    is treated as a cache keyed on the streamwise station x. The `# meta:` line
    records the settings the file was produced with; when they no longer match
    the cache is stale and everything is recomputed, which is what keeps the
    reuse honest.

    Returns ``(cached, x_todo, write_file)``: the rows worth keeping (or None),
    the stations still to compute, and whether the result should be written.

    Cached stations outside ``x_wanted`` are kept rather than dropped, so a run
    over a narrower range never throws away work.
    """
    if cachePolicy not in ("resume", "reuse", "overwrite", "ask"):
        raise ValueError(
            f"cachePolicy must be resume, reuse, overwrite or ask; got {cachePolicy}"
        )

    if cachePolicy == "overwrite" or not os.path.isfile(filename_save):
        return None, x_wanted, True

    if cachePolicy == "ask":
        print(
            f"Max growth rate data already exists at {filename_save}. Overwrite? "
            "(y/N) or skip writing to file (s)"
        )
        choice = input().strip().lower()
        if choice == "y":
            return None, x_wanted, True
        if choice == "s":
            return None, x_wanted, False
        cachePolicy = "reuse"

    try:
        with warnings.catch_warnings():
            # an existing but empty file is a normal outcome of an interrupted
            # run, not something worth warning about
            warnings.simplefilter("ignore")
            cached = np.atleast_2d(np.loadtxt(filename_save))
    except ValueError:
        cached = np.array([])

    if cached.size == 0 or cached.shape[1] < 5:
        print(f"{filename_save} is empty or malformed; recomputing.")
        return None, x_wanted, True

    stored = _parse_meta(filename_save)
    wanted = {k: str(v) for k, v in meta.items()}
    if not stored:
        print(
            f"{filename_save} predates the settings stamp, so it cannot be "
            "checked against the current ones; reusing it anyway."
        )
    else:
        changed = sorted(
            k for k in set(stored) | set(wanted) if stored.get(k) != wanted.get(k)
        )
        if changed:
            print(
                f"Settings changed since {filename_save} was written "
                f"({', '.join(changed)}); recomputing from scratch."
            )
            return None, x_wanted, True

    if cachePolicy == "reuse":
        print(f"Reusing all {len(cached)} cached stations; solver not called.")
        return cached, np.array([], dtype=float), False

    x_cached = cached[:, 4]
    todo = np.array(
        [x for x in x_wanted if not np.any(np.isclose(x_cached, x, atol=X_TOL, rtol=0))],
        dtype=float,
    )

    if todo.size == 0:
        print(f"All {len(x_wanted)} stations already cached; solver not called.")
        return cached, todo, False

    print(f"{len(x_cached)} stations cached, {todo.size} left to compute.")
    return cached, todo, True


def _mergeCache(cached: np.ndarray | None, rows: list[list[float]]) -> np.ndarray:
    """Combine the reused rows with the freshly computed ones, ordered by x."""
    new = np.array(rows, dtype=float).reshape(-1, 5)
    data = new if cached is None else np.vstack((cached, new))

    if data.size == 0:
        raise ValueError(
            "No stations available: nothing was cached and every solver run "
            "failed. Check the swept range and the solver output."
        )

    return data[np.argsort(data[:, 4])]

def saveData(
    filename_save: str,
    reynolds: np.ndarray,
    var_r: np.ndarray,
    ev_r: np.ndarray,
    growthRate: np.ndarray,
    xpositions: np.ndarray,
    temporal: bool,
    meta: dict | None = None,
) -> None:
    var_label = "α" if temporal else "ω"
    ev_label = "ω" if temporal else "α"
    growth_label = "ω_i" if temporal else "-α_i"

    header = (
        f"Re_δ*/1.72 {var_label}_r {ev_label}_r {growth_label} x "
        "# locus of maximum growth rate"
    )
    order = np.argsort(xpositions)
    np.savetxt(
        f"{filename_save}",
        np.column_stack(
            (
                reynolds[order],
                var_r[order],
                ev_r[order],
                growthRate[order],
                xpositions[order],
            )
        ),
        header=header,
        fmt="%f",
    )



def run_analytical_blasius(
    temporal: bool,
    useGaster: bool,
    toml_file: str,
    var_r_min: float,
    var_r_max: float,
    var_r_num: int,
    filename_ev: str,
    re: float,
    gasterOrder: int = GASTER_ORDER,
    gasterLam: float | None = GASTER_SPLINE_LAM,
    cachePolicy: str = CACHE_POLICY,
) -> None:

    extra = "_temporal" if temporal else "_spatial"
    if useGaster:
        extra += f"_gaster{gasterOrder}"

    # reported branch: spatial for the plain spatial run and after Gaster
    reported_temporal = temporal if not useGaster else not temporal
    branch = "temporal" if temporal else "spatial"

    curves_dir = os.path.join(_SCRIPT_DIR, "../../../data/neutralCurves")
    filename_save = os.path.join(
        curves_dir, f"analytical_blasius_maxGrowth{extra}.dat"
    )
    # the neutral curve is order independent (b = 0 wherever omega_i = 0), so it
    # is stored without the Gaster order suffix
    extra_neutral = "_temporal" if temporal else "_spatial"
    if useGaster:
        extra_neutral += "_gaster"
    neutral_curve_file = os.path.join(
        curves_dir, f"analytical_blasius{extra_neutral}.dat"
    )

    x_scan1 = np.arange(-250, -150, 10)
    x_scan2 = np.arange(-150, 1000, 50)
    x_scan3 = np.arange(1000, 3000, 100)
    x_scan = np.sort(np.concatenate((x_scan1, x_scan2, x_scan3)).astype(float))

    meta = {
        "re": re,
        "temporal": temporal,
        "useGaster": useGaster,
        "gasterOrder": gasterOrder,
        "gasterLam": gasterLam,
        "sweep": f"{var_r_min}:{var_r_max}:{var_r_num}",
    }
    cached, x_todo, write_file = resolveCache(
        filename_save, x_scan, meta, cachePolicy
    )
    if not write_file:
        filename_save = ""

    rows: list[list[float]] = []

    if x_todo.size:
        _, ax = plt.subplots()

    for x in x_todo:
        delta_star_ratio = np.sqrt(1.0 + x * BLASIUS_C**2 / re)
        re_solver = re * delta_star_ratio

        print(f"Running for x = {x}, Re_delta* = {re_solver}")

        replacement_line = np.array(
            [
                f"re = {re_solver}",
                'problem = "BoundaryLayer"',
                f'branch = "{branch}"',
                (
                    f"vars_r = {{min = {var_r_min}, "
                    f"max = {var_r_max}, num = {var_r_num}}}"
                ),
            ]
        )
        line_startswith = np.array(["re = ", "problem = ", "branch = ", "vars_r = "])
        editFile(toml_file, replacement_line, line_startswith)

        run_os(toml_file)
        ev_r, ev_i, var_r_solver = read_evs(filename_ev)

        ax.plot(ev_r, ev_i, "o", label=f"x = {x}")

        result = get_maxGrowthRate(
            ev_r, ev_i, var_r_solver, temporal, useGaster, gasterOrder, gasterLam
        )
        if result is None:
            print(f"No usable eigenvalue branch at x = {x}.")
            continue
        var_r_star, ev_r_star, growth_star = result

        rows.append(
            [re_solver / BLASIUS_C, var_r_star, ev_r_star, growth_star, float(x)]
        )

    data = _mergeCache(cached, rows)

    if filename_save:
        saveData(
            filename_save,
            data[:, 0],
            data[:, 1],
            data[:, 2],
            data[:, 3],
            data[:, 4],
            reported_temporal,
            meta,
        )

    plot_max_growthRate(
        data[:, 0],
        data[:, 1],
        data[:, 2],
        data[:, 3],
        data[:, 4],
        reported_temporal,
        neutral_curve_file,
    )


def getXlocation(data_dir: str, pattern: str) -> np.ndarray:
    files = Path(data_dir).glob(pattern)
    x_locations = []
    for file in files:
        name = file.stem
        parts = name.split("_")
        for part in parts:
            if part.startswith("x"):
                try:
                    x_val = float(part[1:])
                    x_locations.append(x_val)
                except ValueError:
                    print(
                        f"Warning: Could not parse x value from '{part}' in filename '{name}'"
                    )
    return np.array(x_locations, dtype=float)


def computeDeltaStar(x: float, re: float, filename: str) -> float:
    analytical_blasius = True

    if analytical_blasius:
        dstar = np.sqrt(1.0 + x * BLASIUS_C**2 / re)
    else:
        # read the velocity profile from the file and compute delta_star
        data = np.loadtxt(filename, skiprows=3)
        y = data[:, 1]
        u = data[:, 2]
        dstar = np.trapezoid(1 - u, y)
    return dstar


def rescaleProfile(filename: str, delta_star: float, saveTo: str) -> None:
    data = np.loadtxt(filename, skiprows=3)
    y = data[:, 1]
    u = data[:, 2]
    y_rescaled = y / delta_star
    data_rescaled = np.column_stack((data[:, 0], y_rescaled, u))
    header = "# x y u"
    np.savetxt(
        saveTo,
        data_rescaled,
        header=header,
        comments="",
        fmt="%f",
    )




def run_sections_baseflow_dns(
    data_dir: str,
    pattern: str,
    temporal: bool,
    useGaster: bool,
    toml_file: str,
    filename_ev: str,
    re: float,
    var_r_min: float,
    var_r_max: float,
    var_r_num: int,
    neutral_curve_file: str | None = None,
    gasterOrder: int = GASTER_ORDER,
    gasterLam: float | None = GASTER_SPLINE_LAM,
    cachePolicy: str = CACHE_POLICY,
) -> None:

    extra = "_temporal" if temporal else "_spatial"
    if useGaster:
        extra += f"_gaster{gasterOrder}"

    # reported branch: spatial for the plain spatial run and after Gaster
    reported_temporal = temporal if not useGaster else not temporal
    var = "Omega" if reported_temporal else "Alpha"
    filename_save = f"{data_dir}/Max{var}_i{extra}.dat"
    branch = "temporal" if temporal else "spatial"

    x_loc = np.sort(getXlocation(data_dir, pattern))
    print(f"Found x locations: {x_loc}")

    _, w = extract_depth_width(data_dir)

    # the three sweep ranges below are part of the result, so they are stamped
    # into the file and a change to any of them invalidates the cache
    sweep_up = (var_r_min, var_r_max, var_r_num)
    sweep_gap = (0.2, 1.5, 60)
    sweep_down = (ALPHA_MIN, ALPHA_MAX, ALPHA_NUM)

    meta = {
        "re": re,
        "temporal": temporal,
        "useGaster": useGaster,
        "gasterOrder": gasterOrder,
        "gasterLam": gasterLam,
        "pattern": pattern,
        "w": w,
        "sweepUp": "{}:{}:{}".format(*sweep_up),
        "sweepGap": "{}:{}:{}".format(*sweep_gap),
        "sweepDown": "{}:{}:{}".format(*sweep_down),
    }
    cached, x_todo, write_file = resolveCache(filename_save, x_loc, meta, cachePolicy)
    if not write_file:
        filename_save = ""

    rows: list[list[float]] = []

    if x_todo.size:
        _, ax = plt.subplots()

    for x in x_todo:
        # substitute * for the actual x in pattern
        file = pattern.replace("*", f"{x}")
        delta_star = computeDeltaStar(x, re, f"{data_dir}/{file}")
        delta_star_ratio = np.sqrt(1.0 + x * BLASIUS_C**2 / re)
        re_local = re * delta_star_ratio
        print(f"Running for x = {x}, Re_delta* = {re_local}")

        re_solver = re_local

        # rescale profile to deltaStar=1
        fileRescaled = "/home/victor/Desktop/orrSommerfeldSolver/data/custom_tmp.dat"
        rescaleProfile(f"{data_dir}/{file}", delta_star, fileRescaled)

        # Written on every station rather than only when crossing into or out of
        # the gap: with stations skipped by the cache, a transition can be missed
        # and the toml would keep whichever range the previous run left behind.
        if 0 < x < w:
            sweep = sweep_gap
        elif x >= w:
            sweep = sweep_down
        else:
            sweep = sweep_up

        replacement_line = np.array(
            [
                f"re = {re_solver}",
                'problem = "Custom"',
                f'branch = "{branch}"',
                f'filenameUprofile = "{fileRescaled}"',
                "vars_r = {{min = {}, max = {}, num = {}}}".format(*sweep),
            ]
        )
        line_startswith = np.array(
            ["re = ", "problem = ", "branch = ", "filenameUprofile = ", "vars_r = "]
        )
        editFile(toml_file, replacement_line, line_startswith)

        run_os(toml_file)
        ev_r, ev_i, var_r_solver = read_evs(filename_ev)

        ax.plot(ev_r, ev_i, "o", label=f"x = {x}")

        result = get_maxGrowthRate(
            ev_r, ev_i, var_r_solver, temporal, useGaster, gasterOrder, gasterLam
        )
        if result is None:
            print(f"No usable eigenvalue branch at x = {x}.")
            continue
        var_r_star, ev_r_star, growth_star = result

        rows.append(
            [re_local / BLASIUS_C, var_r_star, ev_r_star, growth_star, float(x)]
        )

    data = _mergeCache(cached, rows)

    if filename_save:
        saveData(
            filename_save,
            data[:, 0],
            data[:, 1],
            data[:, 2],
            data[:, 3],
            data[:, 4],
            reported_temporal,
            meta,
        )

    plot_max_growthRate(
        data[:, 0],
        data[:, 1],
        data[:, 2],
        data[:, 3],
        data[:, 4],
        reported_temporal,
        neutral_curve_file,
    )


def main() -> None:

    temporal = True
    useGaster = True
    gasterOrder = GASTER_ORDER  # 2 keeps the curvature terms, 1 is classical
    # "resume" only runs the solver for stations not already in the output file;
    # "reuse" never runs it, "overwrite" always does, "ask" prompts
    cachePolicy = CACHE_POLICY

    # var = alpha for temporal stablity analysis and omega for spatial stability analysis
    if temporal:
        ########### for temporal ##############
        var_r_min = 0.01
        var_r_max = 2.5
        var_r_num = 90
    else:
        ########### for spatial ##############
        var_r_min = 0.01
        var_r_max = 0.05
        var_r_num = 30

    assert temporal or not useGaster, (
        "Gaster transformation is only applicable for temporal analysis."
    )

    analytical_blasius = False
    toml_file = "/home/victor/Desktop/orrSommerfeldSolver/config/input.toml"
    filename_ev = "/home/victor/Desktop/orrSommerfeldSolver/data/eigenvalues.dat"
    n = 600
    re = 1000
    setup_toml(toml_file, var_r_min, var_r_max, var_r_num)
    pattern = f"points_n{n}_x*.dat"

    ##############
    # single_case
    ##############

    # case = "d2.25_w31"
    # case = "d1.5_w30"

    # data_dir = "/home/victor/Desktop/PhD/src/bfsRe1000inc/directLinearSolver/blowingSuction/d1.5/data/xSections/"
    # # data_dir = "/home/victor/Desktop/PhD/src/flatPlateRe1000inc/directLinearSolver/blowingSuction/data/"
    # #
    # #
    # if analytical_blasius:
    #     run_analytical_blasius(
    #         temporal,
    #         useGaster,
    #         toml_file,
    #         var_r_min,
    #         var_r_max,
    #         var_r_num,
    #         filename_ev,
    #     )
    #     return
    # else:
    #     run_sections_baseflow_dns(
    #         data_dir, pattern, temporal, useGaster, toml_file, filename_ev
    #     )

    cases = [
        "d0.5_w90",
        "d3_w16",
    ]

    ##############
    # multiple_case
    ##############

    # run_sections_baseflow_dns(
    #     "/home/victor/Desktop/PhD/src/flatPlateRe1000inc/directLinearSolver/blowingSuction/data/",
    #     pattern,
    #     temporal,
    #     useGaster,
    #     toml_file,
    #     filename_ev,
    # )

    if analytical_blasius:
        run_analytical_blasius(
            temporal,
            useGaster,
            toml_file,
            var_r_min,
            var_r_max,
            var_r_num,
            filename_ev,
            re,
            gasterOrder=gasterOrder,
            cachePolicy=cachePolicy,
        )
        return

    for case in cases:
        data_dir = f"/home/victor/Desktop/PhD/src/incGapRe{re}/directLinearSolver/blowingSuction/{case}/data/xSections/"
        # overlaid on the locus when it exists; harmless when it does not
        neutral_curve_file = os.path.join(
            _SCRIPT_DIR,
            f"../../../data/neutralCurves/{case}/neutral_curve_temporal_gaster.dat",
        )
        run_sections_baseflow_dns(
            data_dir,
            pattern,
            temporal,
            useGaster,
            toml_file,
            filename_ev,
            re,
            var_r_min,
            var_r_max,
            var_r_num,
            neutral_curve_file=neutral_curve_file,
            gasterOrder=gasterOrder,
            cachePolicy=cachePolicy,
        )


if __name__ == "__main__":
    main()
