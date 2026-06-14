from pathlib import Path
import os
import subprocess

import matplotlib.pyplot as plt
import numpy as np

from pp.fileManagement import editFile

BLASIUS_C = 1.7207876573
RE_REFERENCE = 1000.0


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


def domega_rdalpha_r(omega_r: np.ndarray, alpha_r: np.ndarray, max_index: int) -> float:
    # Calculate the derivative of omega_r with respect to alpha_r using non-uniform grid (just in case)

    if max_index == 0:
        # Forward difference of order 2
        h = alpha_r[1] - alpha_r[0]
        diff = (-3 * omega_r[0] + 4 * omega_r[1] - omega_r[2]) / (2 * h)
    elif max_index == len(alpha_r) - 1:
        # Backward difference of order 2
        h = alpha_r[-1] - alpha_r[-2]
        diff = (3 * omega_r[-1] - 4 * omega_r[-2] + omega_r[-3]) / (2 * h)
    else:
        h = alpha_r[max_index + 1] - alpha_r[max_index - 1]
        diff = (omega_r[max_index + 1] - omega_r[max_index - 1]) / h

    return diff


def gasterTrans(
    omega_r: np.ndarray, omega_i: np.ndarray, alpha_r: np.ndarray, idx: int
) -> float:
    diff = domega_rdalpha_r(omega_r, alpha_r, idx)
    if diff == 0.0:
        raise ValueError(
            f"d(ω_r)/d(α_r) is zero at index {idx}; cannot apply Gaster transform."
        )
    return -omega_i[idx] / diff


def get_neutral_alpha_r_gaster(
    omega_r: np.ndarray,
    omega_i: np.ndarray,
    alpha_r: np.ndarray,
) -> np.ndarray:
    """
    Estimate the spatial neutral points α_r at α_i = 0 using Gaster transform.

    For each temporal ω_i sign change bracket [j, j+1], compute α_i at both
    endpoints and interpolate α_r at α_i = 0.
    """
    nvals = len(omega_i)
    if nvals < 3:
        return np.array([], dtype=float)

    omega_r_neutral_curve: list[float] = []
    for idx in range(nvals - 1):
        omega_i1 = omega_i[idx]
        omega_i2 = omega_i[idx + 1]

        if omega_i1 * omega_i2 < 0.0:
            alpha_i1 = gasterTrans(omega_r, omega_i, alpha_r, idx)
            alpha_i2 = gasterTrans(omega_r, omega_i, alpha_r, idx + 1)

            if alpha_i1 * alpha_i2 < 0.0:
                omega_r_neutral_curve.append(
                    _interpolate_zero_crossing(
                        omega_r[idx], alpha_i1, omega_r[idx + 1], alpha_i2
                    )
                )

    return np.asarray(omega_r_neutral_curve, dtype=float)


def plot_neutral_curve(
    reynolds: np.ndarray,
    var_r_solver_neutralCurve: np.ndarray,
    xpositions: np.ndarray,
    temporal: bool,
    filename2save: str,
) -> None:
    if reynolds.size == 0 or var_r_solver_neutralCurve.size == 0:
        raise ValueError(
            "No neutral points were found. Check alpha range and solver output."
        )

    _, ax1 = plt.subplots()

    var_label = "α" if temporal else "ω"

    ax1.plot(reynolds, var_r_solver_neutralCurve, "o")
    ax1.set_xlabel("Re_δ*")
    ax1.set_ylabel(var_label + "_r")
    ax1.grid()

    # use the formula x = ((Re_δ*/RE_REFERENCE)^2 - 1) * RE_REFERENCE / BLASIUS_C^2 to get x from Re_δ*
    if xpositions.size > 0 and xpositions.size == reynolds.size:

        def re_to_x(re):
            return (
                ((np.asarray(re) / RE_REFERENCE) ** 2 - 1.0)
                * RE_REFERENCE
                / BLASIUS_C**2
            )

        def x_to_re(x):
            return RE_REFERENCE * np.sqrt(
                1.0 + BLASIUS_C**2 * np.asarray(x) / RE_REFERENCE
            )

        ax2 = ax1.secondary_xaxis(
            "top",
            functions=(re_to_x, x_to_re),
        )

        ax2.set_xlabel("x")

    # save the data for later use
    if filename2save != "":
        np.savetxt(
            f"{filename2save}",
            np.column_stack((reynolds, var_r_solver_neutralCurve, xpositions)),
            header="Re_δ*  "
            + var_label
            + "_r  x # neutral stability curve for "
            + var_label
            + "_i = 0",
            fmt="%f",
        )

    plt.title("Neutral curve " + ("(Temporal)" if temporal else "(Spatial)"))
    plt.tight_layout()
    plt.show()


def run_analytical_blasius(
    temporal: bool,
    gasterTransform: bool,
    toml_file: str,
    var_r_min: float,
    var_r_max: float,
    var_r_num: int,
    filename_ev: str,
) -> None:
    var: list[float] = []
    reynolds: list[float] = []
    xpositions: list[float] = []
    branch = "temporal" if temporal else "spatial"
    filename_save = "../../data/neutral_curve_blasius.dat"

    x_scan = np.arange(-150, 1001, 50)
    var_r_min = var_r_min / 2.0

    for x in x_scan:
        delta_star_ratio = np.sqrt(1.0 + x * BLASIUS_C**2 / RE_REFERENCE)
        re_local = RE_REFERENCE * delta_star_ratio

        # Historical scaling in this workflow:
        # 1) alpha_solver = alpha_ref * delta_star_ratio
        # 2) Re entry in the solver is also scaled by delta_star_ratio
        # 3) with solver output in c, omega_ref = (alpha_solver / delta_star_ratio) * c
        var_scale = delta_star_ratio
        re_solver = re_local * delta_star_ratio

        print(f"Running for x = {x}, Re_delta* = {re_local}")

        replacement_line = np.array(
            [
                f"re = {re_solver}",
                'problem = "BoundaryLayer"',
                f'branch = "{branch}"',
                (
                    f"vars_r = {{min = {var_r_min * var_scale}, "
                    f"max = {var_r_max * var_scale}, num = {var_r_num}}}"
                ),
            ]
        )
        line_startswith = np.array(["re = ", "problem = ", "branch = ", "vars_r = "])
        editFile(toml_file, replacement_line, line_startswith)

        run_os(toml_file)
        ev_r, ev_i, var_r_solver = read_evs(filename_ev)

        # rescale ev_r, ev_i back
        ev_r /= var_scale
        ev_i /= var_scale
        var_r_solver /= var_scale

        if temporal and gasterTransform:
            print(f"Applying Gaster transformation for x = {x}...")
            ev_r_neutral_curve = get_neutral_alpha_r_gaster(ev_r, ev_i, var_r_solver)
        else:
            ev_r_neutral_curve = get_neutral_var_r(ev_r, ev_i)

        if ev_r_neutral_curve.size == 0:
            print(f"No neutral point found at x = {x}.")
            continue

        var.extend(ev_r_neutral_curve.tolist())
        reynolds.extend([re_local] * ev_r_neutral_curve.size)
        xpositions.extend([x] * ev_r_neutral_curve.size)

    plot_neutral_curve(
        np.asarray(reynolds, dtype=float),
        np.asarray(var, dtype=float),
        np.asarray(xpositions, dtype=float),
        temporal if not gasterTransform else not temporal,
        filename_save,
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


def run_sections_baseflow_dns(
    data_dir: str,
    pattern: str,
    temporal: bool,
    gasterTransform: bool,
    toml_file: str,
    filename_ev: str,
) -> None:

    extra = "_temporal" if temporal else "_spatial"
    if gasterTransform:
        extra += "_gaster"
    filename_save = f"{data_dir}/neutral_curve{extra}.dat"
    branch = "temporal" if temporal else "spatial"

    if os.path.isfile(filename_save):
        print(f"Neutral curve data already exists at {filename_save}. Overwrite? (y/N)")
        choice = input().strip().lower()
        if choice != "y":
            print("Plotting existing data without running solver...")
            data = np.loadtxt(filename_save, skiprows=1)
            reynolds = data[:, 0]
            var_r = data[:, 1]
            xpositions = data[:, 2]
            plot_neutral_curve(reynolds, var_r, xpositions, temporal, filename_save)
            return

    x_loc = getXlocation(data_dir, pattern)
    x_loc = np.sort(x_loc)
    print(f"Found x locations: {x_loc}")

    var_r: list[float] = []
    reynolds: list[float] = []
    xpositions: list[float] = []

    for x in x_loc:
        delta_star_ratio = np.sqrt(1.0 + x * BLASIUS_C**2 / RE_REFERENCE)
        re_local = RE_REFERENCE * delta_star_ratio
        print(f"Running for x = {x}, Re_delta* = {re_local}")

        # substitute * for the actual x in pattern
        file = pattern.replace("*", f"{x}")

        replacement_line = np.array(
            [
                f"re = {re_local}",
                'problem = "Custom"',
                f'branch = "{branch}"',
                (f'filenameUprofile = "{data_dir}/{file}"'),
            ]
        )
        line_startswith = np.array(
            ["re = ", "problem = ", "branch = ", "filenameUprofile = "]
        )
        editFile(toml_file, replacement_line, line_startswith)

        run_os(toml_file)
        ev_r, ev_i, var_r_solver = read_evs(filename_ev)

        if temporal and gasterTransform:
            print(f"Applying Gaster transformation for x = {x}...")
            var_r_solver_neutral_curve = get_neutral_alpha_r_gaster(
                ev_r, ev_i, var_r_solver
            )
        else:
            var_r_solver_neutral_curve = get_neutral_var_r(ev_i, var_r_solver)

        if var_r_solver_neutral_curve.size == 0:
            print(f"No neutral point found at x = {x}.")
            continue

        var_r.extend(var_r_solver_neutral_curve.tolist())
        reynolds.extend([re_local] * var_r_solver_neutral_curve.size)
        xpositions.extend([x] * var_r_solver_neutral_curve.size)

    plot_neutral_curve(
        np.asarray(reynolds, dtype=float),
        np.asarray(var_r, dtype=float),
        np.asarray(xpositions, dtype=float),
        temporal if not gasterTransform else not temporal,
        filename_save,
    )


def main() -> None:

    temporal = True
    gasterTransform = True

    # var = alpha for temporal stablity analysis and omega for spatial stability analysis
    if temporal:
        ########### for temporal ##############
        var_r_min = 0.01
        var_r_max = 0.4
        var_r_num = 20
    else:
        ########### for spatial ##############
        var_r_min = 0.01
        var_r_max = 0.2
        var_r_num = 10

    assert temporal or not gasterTransform, (
        "Gaster transformation is only applicable for temporal analysis."
    )

    analytical_blasius = False
    toml_file = "/home/victor/Desktop/orrSommerfeldSolver/config/input.toml"
    filename_ev = "/home/victor/Desktop/orrSommerfeldSolver/data/eigenvalues.dat"
    n = 600

    # case = "d2.25_w31"
    case = "d1.25_w90"
    # data_dir = f"/home/victor/Desktop/PhD/src/incGapRe1000/directLinearSolver/blowingSuctionCoarserMesh/{case}/data/"
    data_dir = "/home/victor/Desktop/PhD/src/flatPlateRe1000inc/directLinearSolver/blowingSuction/wgnInsideDomainDivFree/data/"
    pattern = f"points_n{n}_x*.dat"

    setup_toml(toml_file, var_r_min, var_r_max, var_r_num)

    if analytical_blasius:
        run_analytical_blasius(
            temporal,
            gasterTransform,
            toml_file,
            var_r_min,
            var_r_max,
            var_r_num,
            filename_ev,
        )
    else:
        run_sections_baseflow_dns(
            data_dir, pattern, temporal, gasterTransform, toml_file, filename_ev
        )


if __name__ == "__main__":
    main()
