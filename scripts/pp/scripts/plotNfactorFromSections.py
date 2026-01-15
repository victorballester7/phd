import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.axes import Axes
from pp.codenamesNfactor import code_names
from pp.DeltaN_computation import computeNx, computeDeltaN, interpolate_extrapolate
from pp.colors import colors


def addNxplot(d: float, w: float, x: np.ndarray, Nx: np.ndarray, ax: Axes) -> None:
    """Plot points for Nx curves"""
    ax.plot(x - w, Nx, markersize=2, label=f"d = {d:.2f}, w = {w:.2f}")
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
        "purple": (lambda w: w < 11, "w < 11"),
        "red":    (lambda w: 11 <= w < 21, "11 ≤ w < 21"),
        "blue":   (lambda w: 21 <= w < 31, "21 ≤ w < 31"),
        "orange": (lambda w: 31 <= w < 41, "31 ≤ w < 41"),
        "green":  (lambda w: 41 <= w < 51, "41 ≤ w < 51"),
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

    # script path
    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))

    # output file for deltaN
    filenameDeltaN = "../../../latex/Images/data/DeltaN.dat"
    filenameDeltaN = os.path.join(pathCurrentScript, filenameDeltaN)

    basePath = "../../../src/incNSboeingGapRe1000/directLinearSolver/blowingSuction/"
    basePath = os.path.join(pathCurrentScript, basePath)

    _, ax = plt.subplots(figsize=(8, 6))
    _, ax2 = plt.subplots(figsize=(8, 6))

    # flat plate case
    dataFile_flat = f"../../../src/flatSurfaceRe1000IncNS/directLinearSolver/blowingSuction/wgnInsideDomainDivFree/data/pointsavg_all_n{n}.dat"
    dataFile_flat = os.path.join(pathCurrentScript, dataFile_flat)
    _, _, x_flat, Nx_flat = computeNx(dataFile_flat, doLoo)
    addNxplot(0, 0, x_flat, Nx_flat, ax)

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

    depths = np.array([])
    widths = np.array([])
    deltaNs = np.array([])

    for dw in code_names:
        dataFile_dw = os.path.join(basePath, dw, "data", f"pointsavg_all_n{n}.dat")
        d, w, x, Nx = computeNx(dataFile_dw, doLoo)
        depths = np.append(depths, d)
        widths = np.append(widths, w)
        # compute deltaN
        dN = computeDeltaN(w, x, Nx, x_flat, Nx_flat)
        deltaNs = np.append(deltaNs, dN)


        addNxplot(d, w, x, Nx, ax)
        plotd_vs_dNx(d, w, dN, ax2, plotBFS=False)

    plotd_vs_dNx(0, 0, 0, ax2, plotBFS=True)

    interpolate_extrapolate(depths, widths, deltaNs, filenameDeltaN)

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
