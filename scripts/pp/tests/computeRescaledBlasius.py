import numpy as np

from scripts.plotNeutralCurve import BLASIUS_C
from pp.fileManagement import editFile


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


def main():
    case = "d3_w16"
    x = 7.5
    re = 3000
    filename = f"/home/victor/Desktop/PhD/src/incGapRe{re}/directLinearSolver/blowingSuction/{case}/data/xSections/points_n600_x{x}.dat"
    saveTo = "/home/victor/Desktop/orrSommerfeldSolver/data/custom_tmp.dat"
    filename_toml = "/home/victor/Desktop/orrSommerfeldSolver/config/input.toml"
    dstar = computeDeltaStar(x, re, filename)
    rescaleProfile(filename, dstar, saveTo)
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
            f"re = {re * dstar}",
            "beta = { r = 0.0, i = 0.0 }",
            "useTargetEV = false",
            "vars_r = {min = 0.01, max = 1.15, num = 50}",
            "vars_i = {min = 0.0, max = 0.0, num = 1}",
            'branch = "temporal"',
            "doPlot = true",
            "use_c = false",
            "multipleRun = true",
            "plotUprofile = false",
            "colX = 1",
            "colY = 2",
            "numSkipHeaderLines = 1",
        ]
    )
    editFile(filename_toml, replacement_line, line_startswith)


if __name__ == "__main__":
    main()
