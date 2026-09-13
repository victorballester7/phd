from pathlib import Path
import os
import subprocess

import matplotlib.pyplot as plt
import numpy as np

from pp.fileManagement import editFile, extract_depth_width

BLASIUS_C = 1.7207876573

ALPHA_MIN = 0.01
ALPHA_MAX = 0.5
ALPHA_NUM = 45

OMEGA_MIN = 0.01
OMEGA_MAX = 0.5
OMEGA_NUM = 30


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


def dy_dx(y: np.ndarray, x: np.ndarray, max_index: int) -> float:
    # Calculate the derivative of omega_r with respect to alpha_r using non-uniform grid (just in case)

    if max_index == 0:
        # Forward difference of order 2
        h = x[1] - x[0]
        diff = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h)
    elif max_index == len(x) - 1:
        # Backward difference of order 2
        h = x[-1] - x[-2]
        diff = (3 * y[-1] - 4 * y[-2] + y[-3]) / (2 * h)
    else:
        h = x[max_index + 1] - x[max_index - 1]
        diff = (y[max_index + 1] - y[max_index - 1]) / h

    return diff


def d2y_dx2(y: np.ndarray, x: np.ndarray, max_index: int) -> float:
    # Calculate the second derivative of omega_r with respect to alpha_r using non-uniform grid (just in case)

    if max_index == 0:
        h1 = x[1] - x[0]
        h2 = x[2] - x[1]
        diff2 = ((y[2] - y[1]) / h2 - (y[1] - y[0]) / h1) / ((h1 + h2) / 2)

    elif max_index == len(x) - 1:
        h1 = x[-1] - x[-2]
        h2 = x[-2] - x[-3]
        diff2 = ((y[-1] - y[-2]) / h1 - (y[-2] - y[-3]) / h2) / ((h1 + h2) / 2)
    else:
        h1 = x[max_index] - x[max_index - 1]
        h2 = x[max_index + 1] - x[max_index]

        diff2 = (
            (y[max_index + 1] - y[max_index]) / h2
            - (y[max_index] - y[max_index - 1]) / h1
        ) / ((h1 + h2) / 2)

    return diff2


def gasterTrans(
    omega_r: np.ndarray, omega_i: np.ndarray, alpha_r: np.ndarray, idx: int
) -> tuple[float, float]:
    cg = dy_dx(omega_r, alpha_r, idx)
    order = 1

    if order == 1:
        result = -omega_i[idx] / cg
    else:  # order == 2:
        dcg_dalpha = d2y_dx2(omega_r, alpha_r, idx)
        result = (-omega_i[idx] / cg) * (1 + 0.5 * (dcg_dalpha / cg**2) * omega_i[idx])

    if cg == 0.0:
        raise ValueError(
            f"d(ω_r)/d(α_r) is zero at index {idx}; cannot apply Gaster transform."
        )
    return -omega_i[idx] / cg, cg


def get_neutral_alpha_r_gaster(
    omega_r: np.ndarray,
    omega_i: np.ndarray,
    alpha_r: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Estimate the spatial neutral points α_r at α_i = 0 using Gaster transform.

    For each temporal ω_i sign change bracket [j, j+1], compute α_i at both
    endpoints and interpolate α_r at α_i = 0.
    """
    # nvals = len(omega_i)
    # if nvals < 3:
    #     return np.array([], dtype=float), np.array([], dtype=float)
    #
    # omega_r_neutral_curve: list[float] = []
    # cg_neutral: list[float] = []
    # for idx in range(nvals - 1):
    #     omega_i1 = omega_i[idx]
    #     omega_i2 = omega_i[idx + 1]
    #
    #     if omega_i1 * omega_i2 < 0.0:
    #         alpha_i1, cg1 = gasterTrans(omega_r, omega_i, alpha_r, idx)
    #         alpha_i2, cg2 = gasterTrans(omega_r, omega_i, alpha_r, idx + 1)
    #
    #         if alpha_i1 * alpha_i2 < 0.0:
    #             omega_r_neutral_curve.append(
    #                 _interpolate_zero_crossing(
    #                     omega_r[idx], alpha_i1, omega_r[idx + 1], alpha_i2
    #                 )
    #             )
    #             cg_neutral.append(
    #                 _interpolate_zero_crossing(cg1, alpha_i1, cg2, alpha_i2)
    #             )
    #
    # return np.asarray(omega_r_neutral_curve, dtype=float), np.asarray(
    #     cg_neutral, dtype=float
    # )
    #

    cg = np.gradient(omega_r, alpha_r)
    cg = np.abs(cg)  # take absolute value to avoid negative group velocity (and avoid numerical misconvergence)
    alpha_i = -omega_i / cg
    neutral_indices = np.where(np.diff(np.sign(alpha_i)))[0]
    omega_r_neutral_curve: np.ndarray = np.array([], dtype=float)
    alpha_r_neutral_curve: np.ndarray = np.array([], dtype=float)
    cg_neutral: np.ndarray = np.array([], dtype=float)

    for idx in neutral_indices:
        alpha_i1 = alpha_i[idx]
        alpha_i2 = alpha_i[idx + 1]
        if alpha_i1 * alpha_i2 < 0.0:
            omega_r_neutral_curve = np.append(
                omega_r_neutral_curve,
                _interpolate_zero_crossing(
                    omega_r[idx], alpha_i1, omega_r[idx + 1], alpha_i2
                ),
            )
            alpha_r_neutral_curve = np.append(
                alpha_r_neutral_curve,
                _interpolate_zero_crossing(
                    alpha_r[idx], alpha_i1, alpha_r[idx + 1], alpha_i2
                ),
            )
            cg_neutral = np.append(
                cg_neutral,
                _interpolate_zero_crossing(cg[idx], alpha_i1, cg[idx + 1], alpha_i2),
            )

    return omega_r_neutral_curve, alpha_r_neutral_curve, cg_neutral


def plot_neutral_curve(
    reynolds: np.ndarray,
    var_r_solver_neutralCurve: np.ndarray,
    ev_r_solver_neutralCurve: np.ndarray,
    xpositions: np.ndarray,
    temporal: bool,
) -> None:
    if reynolds.size == 0 or var_r_solver_neutralCurve.size == 0:
        raise ValueError(
            "No neutral points were found. Check alpha range and solver output."
        )

    _, ax1 = plt.subplots()
    _, ax2 = plt.subplots()

    var_label = "α" if temporal else "ω"
    ev_label = "ω" if temporal else "α"
    reynolds = 1.72 * reynolds  # rescale to Re_δ* = 1.72 * Re_δ*
    ax1.plot(reynolds, var_r_solver_neutralCurve, "o")
    ax2.plot(reynolds, ev_r_solver_neutralCurve, "o")
    ax1.set_xlabel("Re_δ*")
    ax2.set_xlabel("Re_δ*")
    ax1.set_ylabel(var_label + "_r")
    ax2.set_ylabel(ev_label + "_r")
    ax1.set_title("Neutral curve " + ("(Temporal)" if temporal else "(Spatial)"))
    ax2.set_title("Neutral curve " + ("(Temporal)" if temporal else "(Spatial)"))
    ax1.grid()
    ax2.grid()

    plt.tight_layout()
    plt.show()


def saveData(
    filename_save: str,
    reynolds: np.ndarray,
    var_r_solver_neutralCurve: np.ndarray,
    ev_r_solver_neutralCurve: np.ndarray,
    xpositions: np.ndarray,
    temporal: bool,
) -> None:
    var_label = "α" if temporal else "ω"
    ev_label = "ω" if temporal else "α"
    np.savetxt(
        f"{filename_save}",
        np.column_stack(
            (
                reynolds,
                var_r_solver_neutralCurve,
                ev_r_solver_neutralCurve,
                xpositions,
            )
        ),
        header="Re_δ*/1.72 "
        + var_label
        + "_r "
        + ev_label
        + "_r x # neutral stability curve for "
        + var_label
        + "_i = 0",
        fmt="%f",
    )


def run_analytical_blasius(
    temporal: bool,
    gasterTransform: bool,
    toml_file: str,
    var_r_min: float,
    var_r_max: float,
    var_r_num: int,
    filename_ev: str,
    re: float,
) -> None:

    var_r: np.ndarray = np.array([], dtype=float)
    evs_r: np.ndarray = np.array([], dtype=float)
    reynolds: np.ndarray = np.array([], dtype=float)
    xpositions: np.ndarray = np.array([], dtype=float)

    extra = "_temporal" if temporal else "_spatial"
    if gasterTransform:
        extra += "_gaster"
    filename_save = f"../../data/neutralCurves/analytical_blasius{extra}.dat"
    branch = "temporal" if temporal else "spatial"

    if os.path.isfile(filename_save):
        print(
            f"Neutral curve data already exists at {filename_save}. Overwrite? (y/N) or skip writon to File (s)"
        )
        choice = input().strip().lower()
        if choice != "y" and choice != "s":
            print("Plotting existing data without running solver...")
            data = np.loadtxt(filename_save, skiprows=1)
            reynolds = data[:, 0]
            var_r = data[:, 1]
            evs_r = data[:, 2]
            xpositions = data[:, 3]
            plot_neutral_curve(reynolds, var_r, evs_r, xpositions, temporal if not gasterTransform else not temporal)
            return
        if choice == "s":
            filename_save = ""

    x_scan1 = np.arange(-250, -150, 10)
    x_scan2 = np.arange(-150, 1000, 50)
    x_scan3 = np.arange(1000, 3000, 100)
    x_scan = np.concatenate((x_scan1, x_scan2, x_scan3))

    _, ax = plt.subplots()
    _, ax2 = plt.subplots()
    for x in x_scan:
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

        if temporal and gasterTransform:
            print(f"Applying Gaster transformation for x = {x}...")
            var_r_solver_neutral_curve, ev_r_solver_neutral_curve, cg_neutral = (
                get_neutral_alpha_r_gaster(ev_r, ev_i, var_r_solver)
            )
            evs_r = np.append(evs_r, ev_r_solver_neutral_curve)
            # cg.append(cg_neutral.tolist())
            # ax2.plot([re_solver] * len(cg_neutral), cg_neutral, "o", label=f"x = {x}")
        else:
            var_r_solver_neutral_curve = get_neutral_var_r(ev_i, var_r_solver)

        if var_r_solver_neutral_curve.size == 0:
            print(f"No neutral point found at x = {x}.")
            continue

        var_r = np.append(var_r, var_r_solver_neutral_curve)
        reynolds = np.append(reynolds, [re_solver] * var_r_solver_neutral_curve.size)
        xpositions = np.append(xpositions, [x] * var_r_solver_neutral_curve.size)
    ax.legend()

    reynolds /= BLASIUS_C

    if filename_save:
        saveData(
            filename_save,
            np.asarray(reynolds, dtype=float),
            np.asarray(var_r, dtype=float),
            np.asarray(evs_r, dtype=float),
            np.asarray(xpositions, dtype=float),
            temporal if not gasterTransform else not temporal,
        )

    plot_neutral_curve(
        np.asarray(reynolds, dtype=float),
        np.asarray(var_r, dtype=float),
        np.asarray(evs_r, dtype=float),
        np.asarray(xpositions, dtype=float),
        temporal if not gasterTransform else not temporal,
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


def computeDeltaStar(x: float, re:float, filename: str) -> float:
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
    gasterTransform: bool,
    toml_file: str,
    filename_ev: str,
    re: float,
) -> None:

    extra = "_temporal" if temporal else "_spatial"
    if gasterTransform:
        extra += "_gaster"
    filename_save = f"{data_dir}/neutral_curve{extra}.dat"
    branch = "temporal" if temporal else "spatial"

    if os.path.isfile(filename_save):
        print(
            f"Neutral curve data already exists at {filename_save}. Overwrite? (y/N) or skip writon to File (s)"
        )
        choice = input().strip().lower()
        if choice != "y" and choice != "s":
            print("Plotting existing data without running solver...")
            data = np.loadtxt(filename_save, skiprows=1)
            reynolds = data[:, 0]
            var_r = data[:, 1]
            ev_r = data[:, 2]
            xpositions = data[:, 3]
            plot_neutral_curve(reynolds, var_r, ev_r, xpositions, temporal)
            return
        if choice == "s":
            filename_save = ""

    x_loc = getXlocation(data_dir, pattern)
    x_loc = np.sort(x_loc)
    print(f"Found x locations: {x_loc}")

    var_r: np.ndarray = np.array([], dtype=float)
    evs_r: np.ndarray = np.array([], dtype=float)
    reynolds: np.ndarray = np.array([], dtype=float)
    xpositions: np.ndarray = np.array([], dtype=float)

    _, ax = plt.subplots()

    _, w = extract_depth_width(data_dir)
    
    x_in_gap_previous = False
    for x in x_loc:
        # if x < 0 or x > 30:
        #     print(f"Skipping x = {x}")
        #     continue

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

        replacement_line = np.array(
            [
                f"re = {re_solver}",
                'problem = "Custom"',
                f'branch = "{branch}"',
                f'filenameUprofile = "{fileRescaled}"',
                # (f'filenameUprofile = "{data_dir}/{file}"'),
            ]
        )
        line_startswith = np.array(
            ["re = ", "problem = ", "branch = ", "filenameUprofile = "]
        )
        if x > 0 and x < w:
            x_in_gap_previous = True
            replacement_line = np.append(
                replacement_line,
                "vars_r = {min = 0.2, max = 1.5, num = 60}",
            )
            line_startswith = np.append(line_startswith, "vars_r = ")

        if x_in_gap_previous and x > w:
            x_in_gap_previous = False
            replacement_line = np.append(
                replacement_line,
                f"vars_r = {{min = {ALPHA_MIN}, max = {ALPHA_MAX}, num = {ALPHA_NUM}}}",
            )
            line_startswith = np.append(line_startswith, "vars_r = ")

        editFile(toml_file, replacement_line, line_startswith)

        run_os(toml_file)
        ev_r, ev_i, var_r_solver = read_evs(filename_ev)

        # rescale back the eigenvalues to the original scale for plotting
        # var_r_solver /= delta_star_ratio
        # ev_r /= delta_star_ratio
        # ev_i /= delta_star_ratio

        ax.plot(ev_r, ev_i, "o", label=f"x = {x}")

        if temporal and gasterTransform:
            print(f"Applying Gaster transformation for x = {x}...")
            var_r_solver_neutral_curve, ev_r_solver_neutral_curve, cg_neutral = (
                get_neutral_alpha_r_gaster(ev_r, ev_i, var_r_solver)
            )
            evs_r = np.append(evs_r, ev_r_solver_neutral_curve)
        else:
            var_r_solver_neutral_curve = get_neutral_var_r(ev_i, var_r_solver)

        if var_r_solver_neutral_curve.size == 0:
            print(f"No neutral point found at x = {x}.")
            continue

        var_r = np.append(var_r, var_r_solver_neutral_curve)
        reynolds = np.append(reynolds, [re_local] * var_r_solver_neutral_curve.size)
        xpositions = np.append(xpositions, [x] * var_r_solver_neutral_curve.size)
    ax.legend()
    reynolds /= BLASIUS_C

    if filename_save:
        saveData(
            filename_save,
            reynolds,
            var_r,
            evs_r,
            xpositions,
            temporal if not gasterTransform else not temporal,
        )

    plot_neutral_curve(
        reynolds,
        var_r,
        evs_r,
        xpositions,
        temporal if not gasterTransform else not temporal,
    )


def main() -> None:

    temporal = True
    gasterTransform = True

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

    assert temporal or not gasterTransform, (
        "Gaster transformation is only applicable for temporal analysis."
    )

    analytical_blasius = True
    toml_file = "/home/victor/Desktop/orrSommerfeldSolver/config/input.toml"
    filename_ev = "/home/victor/Desktop/orrSommerfeldSolver/data/eigenvalues.dat"
    n = 600
    re = 3000
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
    #         gasterTransform,
    #         toml_file,
    #         var_r_min,
    #         var_r_max,
    #         var_r_num,
    #         filename_ev,
    #     )
    #     return
    # else:
    #     run_sections_baseflow_dns(
    #         data_dir, pattern, temporal, gasterTransform, toml_file, filename_ev
    #     )


    cases = [
        # "d1.5_w10",
        # "d1.5_w15",
        # "d1.5_w20",
        # "d1.5_w25",
        # "d1.5_w30",
        # "d1.5_w35",
        # "d1.5_w40",
        # "d1.5_w45",
        # "d1.5_w70",
        # "d1.75_w33",
        # "d2_w41",
        # "d3.5_w21",
        "d1.25_w38",
        "d0.75_w93",
        "d3_w14",
        "d1.5_w33",
        # "d0.25_w40",
        # "d0.75_w40",
        # "d1_w39",
        # "d1.25_w40",
        # "d1.5_w40",
        # "d1.75_w40",
        # "d1_w20",
        # "d2_w20",
        # "d2.5_w20",
        # "d3_w21",
    ]

    ##############
    # multiple_case
    ##############

    # run_sections_baseflow_dns(
    #     "/home/victor/Desktop/PhD/src/flatPlateRe1000inc/directLinearSolver/blowingSuction/data/",
    #     pattern,
    #     temporal,
    #     gasterTransform,
    #     toml_file,
    #     filename_ev,
    # )

    for case in cases:
        data_dir = f"/home/victor/Desktop/PhD/src/incGapRe3000/directLinearSolver/blowingSuction/{case}/data/xSections/"
        run_sections_baseflow_dns(
            data_dir, pattern, temporal, gasterTransform, toml_file, filename_ev, re
        )


if __name__ == "__main__":
    main()
