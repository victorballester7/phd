import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.axes import Axes
from pp.codenamesNfactor import code_names
from pp.DeltaN_computation import (
    computeNx,
    computeDeltaN,
    gpr_interpolate,
    write_deltaN_file,
    plot_deltaN_grid,
)
from pp.colors import colors
from pp.fileManagement import extract_width_depth

DEPTHS = np.array(
    [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4]
)
LINES = [
    "-",
    "--",
    "-.",
    ":",
    (0, (1, 1)),
    (0, (5, 1)),
    (0, (3, 1, 1, 1)),
    (0, (3, 1, 1, 1, 1, 1)),
    (5, (10, 3)),
]


def addNxplot(
    d: float, w: float, x: np.ndarray, Nx: np.ndarray, ax: Axes, c, ls, linewidth=1
) -> None:
    """Plot points for Nx curves"""
    ax.plot(
        x - w,
        Nx,
        markersize=2,
        label=f"d = {d:.2f}, w = {w:.2f}",
        color=c,
        linestyle=ls,
    )
    return


def plotd_vs_dNx(d: float, w: float, dN: float, ax: Axes, plotBFS: bool) -> None:
    """Plot d vs DeltaN(x)"""

    if not hasattr(plotd_vs_dNx, "used_labels"):
        plotd_vs_dNx.used_labels = set()

    if plotBFS:
        xbfs = np.linspace(0, 1.25, 100)
        ybfs = 4.4 * xbfs
        ax.plot(xbfs, ybfs, "--", color="black", label="BFS (Crouch et al. 2006)")
        return

    color_map = {
        "tab:blue": (lambda w: w < 11, "w < 11"),
        "tab:orange": (lambda w: 11 <= w < 21, "11 ≤ w < 21"),
        "tab:green": (lambda w: 21 <= w < 31, "21 ≤ w < 31"),
        "tab:red": (lambda w: 31 <= w < 41, "31 ≤ w < 41"),
        "tab:purple": (lambda w: 41 <= w < 51, "41 ≤ w < 51"),
        "tab:brown": (lambda w: 51 <= w < 61, "51 ≤ w < 61"),
        "tab:pink": (lambda w: 61 <= w, "w ≥ 61"),
    }

    col, label = None, None
    for c, (cond, lbl) in color_map.items():
        if cond(w):
            col, label = c, lbl
            break

    if col is None:
        print(colors.WARNING + f"Width {w} too large, skipping" + colors.ENDC)
        return

    if label not in plotd_vs_dNx.used_labels:
        ax.plot(d, dN, "o", color=col, label=label)
        plotd_vs_dNx.used_labels.add(label)
    else:
        ax.plot(d, dN, "o", color=col)


def main():
    """
    Read the data files that contain the time averaged reynolds stresses for several x and y locations, computes the amplitude in the y direction (L2 or Loo) and then computes the N factor as Nx = log(A/A0) where A0 is the amplitude at the first x location.
    Plots the N factor for several cases (gap configurations) in the same plot as well as the flat plate case for comparison.
    """

    # Parameters to change
    n = 600  # number of points in the wall normal direction
    doLoo = False  # if False, L2 norm is computed
    plotBFSdata = False
    plotdeepGapdata = True

    assert not (plotBFSdata and plotdeepGapdata), "Cannot plot BFS data and deep gap data at the same time, as they have different line styles and colors for the same depths, which makes the plot confusing. Please choose one of the two options."

    # script path
    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))

    # output file for deltaN
    filenameDeltaN = "../../../latex/Images/data/DeltaN.dat"
    filenameDeltaN = os.path.join(pathCurrentScript, filenameDeltaN)

    basePath = "../../../src/incGapRe1000/directLinearSolver/blowingSuction/"
    basePath = os.path.join(pathCurrentScript, basePath)

    _, ax = plt.subplots(figsize=(8, 6))
    _, ax2 = plt.subplots(figsize=(8, 6))

    # flat plate case
    dataFile_flat = f"../../../src/flatPlateRe1000inc/directLinearSolver/blowingSuction/wgnInsideDomainDivFree/data/pointsavg_n{n}.dat"
    dataFile_flat = os.path.join(pathCurrentScript, dataFile_flat)
    x_flat, Nx_flat = computeNx(dataFile_flat, doLoo)
    addNxplot(0, 0, x_flat, Nx_flat, ax, "black", "-", linewidth=2)

    # just to do a quick check
    # code_names2 = [
    #     "d2_w10",
    #     "d2_w102",
    #     "d2_w20",
    #     "d2_w202",
    #     "d2_w32",
    #     "d2_w322",
    #     "d3.25_w10",
    #     "d3.25_w102",
    # ]
    
    basePathBFS = "../../../src/bfsRe1000inc/directLinearSolver/blowingSuction/"
    basePathBFS = os.path.join(pathCurrentScript, basePathBFS)
    bfs_codenames = [
        "d0.25",
        "d0.5",
        "d0.75",
        "d1",
        "d1.25",
    ]

    basePathDeepGap = "../../../src/deepGapRe1000inc/directLinearSolver/blowingSuction/"
    basePathDeepGap = os.path.join(pathCurrentScript, basePathDeepGap)
    deepGap_codenames = [
        "d100_w5",
        "d100_w10",
        "d100_w15"
    ]

    if plotBFSdata:
        for dw in bfs_codenames:
            dataFile_dw = os.path.join(basePathBFS, dw, "data", f"pointsavg_n{n}.dat")
            d, w = extract_width_depth(dataFile_dw)
            x, Nx = computeNx(dataFile_dw, doLoo)
            # get the color in tab20 based on the index of the depth in DEPTHS
            idx = np.argmin(np.abs(DEPTHS - d))
            c = plt.get_cmap("tab20")(idx)
            ls = LINES[idx % len(LINES)]
            addNxplot(d, w, x, Nx, ax, c, ls)
            dN = computeDeltaN(w, x, Nx, x_flat, Nx_flat)
            plotd_vs_dNx(d, w, dN, ax2, plotBFS=False)

    if plotdeepGapdata:
        ls_styles_deepGap = [":", "--", "-."]  # line styles for deep gap cases
        for dw in deepGap_codenames:
            dataFile_dw = os.path.join(basePathDeepGap, dw, "data", f"pointsavg_n{n}.dat")
            d, w = extract_width_depth(dataFile_dw)
            x, Nx = computeNx(dataFile_dw, doLoo)
            # get the color in tab20 based on the index of the depth in DEPTHS
            idx = np.argmin(np.abs(DEPTHS - d))
            c = "black"
            ls = ls_styles_deepGap[deepGap_codenames.index(dw) % len(ls_styles_deepGap)]  # line style based on the index of the deep gap case
            addNxplot(d, w, x, Nx, ax, c, ls)
            dN = computeDeltaN(w, x, Nx, x_flat, Nx_flat)
            plotd_vs_dNx(d, w, dN, ax2, plotBFS=False)

    depths = np.array([])
    widths = np.array([])
    deltaNs = np.array([])

    prevD = None
    idx_ls = 0
    for dw in code_names:
        dataFile_dw = os.path.join(basePath, dw, "data", f"pointsavg_n{n}.dat")
        d, w = extract_width_depth(dataFile_dw)

        ############## TEMPORARY ##############
        if plotBFSdata and d > 1.25:
            # skip large depths for now, as they are not relevant for the report and take a long time to compute
            print(colors.WARNING + f"Depth {d} too large, skipping" + colors.ENDC)
            continue

        if plotdeepGapdata and w > 20:
            # skip large widths for now, as they are not relevant for the report and take a long time to compute
            print(colors.WARNING + f"Width {w} too large, skipping" + colors.ENDC)
            continue

        #######################################

        x, Nx = computeNx(dataFile_dw, doLoo)
        depths = np.append(depths, d)
        widths = np.append(widths, w)
        # compute deltaN
        dN = computeDeltaN(w, x, Nx, x_flat, Nx_flat)
        deltaNs = np.append(deltaNs, dN)

        # get the color in tab20 based on the index of the depth in DEPTHS
        idx = np.argmin(np.abs(DEPTHS - d))
        c = plt.get_cmap("tab20")(idx)

        if prevD is not None and d == prevD:
            idx_ls += 1
        else:
            idx_ls = 0
            prevD = d

        idx_tmp = (idx_ls % len(LINES)) if not plotBFSdata else ((idx_ls % len(LINES) - 1) + 1)  # skip the first line style for the first case of each depth, as it is already used for the BFS case
        ls = LINES[idx_tmp] # skip the first line style for the first case of each depth, as it is already used for the BFS case
        addNxplot(d, w, x, Nx, ax, c, ls)
        plotd_vs_dNx(d, w, dN, ax2, plotBFS=False)

    plotd_vs_dNx(0, 0, 0, ax2, plotBFS=True)

    # GPR interpolation and plotting (separated)
    gpr_result = gpr_interpolate(depths, widths, deltaNs)
    write_deltaN_file(gpr_result, filenameDeltaN)
    plot_deltaN_grid(gpr_result)

    ax.set_title(f"n(x) factor using {'Loo' if doLoo else 'L2'} norm")
    ax.set_xlabel("x - w")
    ax.set_ylabel("n(x)")
    ax.grid()
    ax.legend()
    ax2.set_title("Delta N as function of d for different w")
    ax2.set_xlabel("d")
    ax2.set_ylabel("Delta N")
    ax2.grid()
    ax2.legend()

    plt.show()


if __name__ == "__main__":
    main()
