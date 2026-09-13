import numpy as np
import matplotlib.pyplot as plt
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from matplotlib.axes import Axes
from pp.codenamesNfactor import code_names_Re
from pp.DeltaN_computation import (
    computeNx,
    computeDeltaN,
    gpr_interpolateWithRe,
    write_deltaN_file,
    plot_deltaN_grid,
)
from pp.colors import colors
from pp.fileManagement import extract_depth_width

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

_worker_x_flat = None
_worker_Nx_flat = None


def init_nfactor_worker(x_flat, Nx_flat):
    global _worker_x_flat, _worker_Nx_flat
    _worker_x_flat = x_flat
    _worker_Nx_flat = Nx_flat


def compute_nfactor_case(data_file, doLoo):
    """Load one gap case and compute its N factor and Delta N independently."""
    d, w = extract_depth_width(data_file)
    x, Nx = computeNx(data_file, doLoo)
    if len(x) == 0 or len(Nx) == 0:
        return d, w, x, Nx, np.nan

    dN = computeDeltaN(w, x, Nx, _worker_x_flat, _worker_Nx_flat)
    return d, w, x, Nx, dN


def addNxplot(
    d: float, w: float, x: np.ndarray, Nx: np.ndarray, ax: Axes, c, ls, linewidth=1
) -> None:
    """Plot points for Nx curves"""

    ax.plot(
        x,
        Nx,
        markersize=2,
        label=f"d = {d:.2f}, w = {w:.2f}",
        color=c,
        linestyle=ls,
        linewidth=linewidth,
    )
    

    x_start = w + 80
    x_end = x_start + 75
    idx_start = np.argmin(np.abs(x - x_start))
    idx_end = np.argmin(np.abs(x - x_end))
    ax.plot([x_start], [Nx[idx_start]], "o", color=c, markersize=5)
    ax.plot([x_end], [Nx[idx_end]], "x", color=c, markersize=5)

    return


def plotd_vs_dNx(
    d: float, w: float, dN: float, ax: Axes, plotBFS: bool, marker: str = "o"
) -> None:
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

    if marker == "s":
        label = "BFS"
    if label not in plotd_vs_dNx.used_labels:
        ax.plot(d, dN, marker, color=col, label=label)
        plotd_vs_dNx.used_labels.add(label)
    else:
        ax.plot(d, dN, marker, color=col)


def main():
    """
    Read the data files that contain the time averaged reynolds stresses for several x and y locations, computes the amplitude in the y direction (L2 or Loo) and then computes the N factor as Nx = log(A/A0) where A0 is the amplitude at the first x location.
    Plots the N factor for several cases (gap configurations) in the same plot as well as the flat plate case for comparison.
    """

    # Parameters to change
    n = 600  # number of points in the wall normal direction
    doLoo = False  # if False, L2 norm is computed
    re = 800  # Reynolds slice to predict/write from the multi-Re GPR model
    plotBFSandFFSdata = False
    plotdeepGapdata = False
    useOnly1Re = True

    if re != 1000:
        plotBFSandFFSdata = False
        plotdeepGapdata = False

    assert not (plotBFSandFFSdata and plotdeepGapdata), (
        "Cannot plot BFS data and deep gap data at the same time, as they have different line styles and colors for the same depths, which makes the plot confusing. Please choose one of the two options."
    )

    # script path
    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))

    # output file for deltaN
    filenameDeltaN = f"../../../data/deltaN/DeltaN{re}.dat"
    filenameDeltaN = os.path.join(pathCurrentScript, filenameDeltaN)

    _, ax = plt.subplots(figsize=(8, 6))
    _, ax2 = plt.subplots(figsize=(8, 6))

    def flat_plate_file(reynolds: int) -> str:
        relative_path =  f"../../../src/flatPlateRe{reynolds}inc/directLinearSolver/blowingSuction/data/pointsavg_n{n}.dat"
        return os.path.join(pathCurrentScript, relative_path)

    # Flat plate curve for the target slice plot.
    dataFile_flat = flat_plate_file(re)
    x_flat, Nx_flat = computeNx(dataFile_flat, doLoo)
    # print(x_flat, Nx_flat)
    addNxplot(0, 0, x_flat, Nx_flat, ax, "black", "-", linewidth=3)
    flat_plate_cache = {re: (x_flat, Nx_flat)}

    basePathBFS = f"../../../src/bfsRe{re}inc/directLinearSolver/blowingSuction/"
    basePathBFS = os.path.join(pathCurrentScript, basePathBFS)
    bfs_codenames = [
        "d0.25",
        "d0.5",
        "d0.75",
        "d1",
        "d1.25",
        "d1.5",
    ]

    basePathFFS = "../../../src/ffsRe1000inc/directLinearSolver/blowingSuction/"
    basePathFFS = os.path.join(pathCurrentScript, basePathFFS)
    ffs_codenames = [
        "d0.5",
        "d1",
        "d1.5",
    ]

    basePathDeepGap = "../../../src/deepGapRe1000inc/directLinearSolver/blowingSuction/"
    basePathDeepGap = os.path.join(pathCurrentScript, basePathDeepGap)
    deepGap_codenames = ["d100_w5", "d100_w10", "d100_w15"]

    bfs_data = []

    if plotBFSandFFSdata:
        for dw in bfs_codenames:
            dataFile_dw = os.path.join(basePathBFS, dw, "data", f"pointsavg_n{n}.dat")
            d, w = extract_depth_width(dataFile_dw)
            x, Nx = computeNx(dataFile_dw, doLoo)
            # get the color in tab20 based on the index of the depth in DEPTHS
            idx = np.argmin(np.abs(DEPTHS - d))
            c = plt.get_cmap("tab20")(idx)
            ls = LINES[idx % len(LINES)]
            addNxplot(d, w, x, Nx, ax, c, ls, linewidth=3)
            dN = computeDeltaN(w, x, Nx, x_flat, Nx_flat, x_start=150, x_end=350)
            bfs_data.append([d, dN])
            plotd_vs_dNx(d, w, dN, ax2, plotBFS=False, marker="s")

        # for dw in ffs_codenames:
        #     dataFile_dw = os.path.join(basePathFFS, dw, "data", f"pointsavg_n{n}.dat")
        #     d, w = extract_width_depth(dataFile_dw)
        #     x, Nx = computeNx(dataFile_dw, doLoo)
        #     # get the color in tab20 based on the index of the depth in DEPTHS
        #     idx = np.argmin(np.abs(DEPTHS - d))
        #     c = plt.get_cmap("tab20")(idx)
        #     ls = LINES[idx % len(LINES)]
        #     addNxplot(d, w, x, Nx, ax, c, ls, linewidth=3)
        #     dN = computeDeltaN(w, x, Nx, x_flat, Nx_flat, x_start=150, x_end=350)
        #     ffs_data.append([d,dN])
        #     plotd_vs_dNx(d, w, dN, ax2, plotBFS=False)

        bfs_data = np.array(bfs_data)

        # fit a line to the bfs data
        xx = np.linspace(0.18, 1.3, 100)
        coeffs = np.polyfit(bfs_data[:, 0], bfs_data[:, 1], 1)
        yy = np.polyval(coeffs, xx)
        ax2.plot(xx, yy, "--", color="purple", label="BFS fit")

    if plotdeepGapdata:
        ls_styles_deepGap = [":", "--", "-."]  # line styles for deep gap cases
        for dw in deepGap_codenames:
            dataFile_dw = os.path.join(
                basePathDeepGap, dw, "data", f"pointsavg_n{n}.dat"
            )
            d, w = extract_depth_width(dataFile_dw)
            x, Nx = computeNx(dataFile_dw, doLoo)
            # get the color in tab20 based on the index of the depth in DEPTHS
            idx = np.argmin(np.abs(DEPTHS - d))
            c = "black"
            ls = ls_styles_deepGap[
                deepGap_codenames.index(dw) % len(ls_styles_deepGap)
            ]  # line style based on the index of the deep gap case
            addNxplot(d, w, x, Nx, ax, c, ls, linewidth=2)
            dN = computeDeltaN(w, x, Nx, x_flat, Nx_flat)
            plotd_vs_dNx(d, w, dN, ax2, plotBFS=False)

    depths = []
    widths = []
    reynolds = []
    deltaNs = []

    prevD = None
    idx_ls = 0
    for r in code_names_Re.keys():
        if useOnly1Re and r != re:
            continue
        if r not in flat_plate_cache:
            flat_plate_cache[r] = computeNx(flat_plate_file(r), doLoo)

        x_flat_r, Nx_flat_r = flat_plate_cache[r]
        if len(x_flat_r) == 0 or len(Nx_flat_r) == 0:
            print(
                colors.WARNING
                + f"Skipping Re={r} because the flat plate reference could not be read"
                + colors.ENDC
            )
            continue

        basePath = f"../../../src/incGapRe{r}/directLinearSolver/blowingSuction/"
        basePath = os.path.join(pathCurrentScript, basePath)

        data_files = [
            os.path.join(basePath, dw, "data", f"pointsavg_n{n}.dat")
            for dw in code_names_Re[r]
        ]
        eligible_files = []
        for data_file in data_files:
            _, w = extract_depth_width(data_file)
            if plotdeepGapdata and w > 20:
                print(colors.WARNING + f"Width {w} too large, skipping" + colors.ENDC)
                continue
            eligible_files.append(data_file)

        if not eligible_files:
            continue

        with ProcessPoolExecutor(
            max_workers=min(len(eligible_files), os.cpu_count() or 1),
            initializer=init_nfactor_worker,
            initargs=(x_flat_r, Nx_flat_r),
        ) as executor:
            results = executor.map(
                compute_nfactor_case,
                eligible_files,
                repeat(doLoo),
                chunksize=1,
            )

        for d, w, x, Nx, dN in results:
            if len(x) == 0 or len(Nx) == 0:
                print(
                    colors.WARNING
                    + "Skipping case due to error in reading data"
                    + colors.ENDC
                )
                continue
            depths.append(d)
            widths.append(w)
            reynolds.append(r)
            deltaNs.append(dN)

            # get the color in tab20 based on the index of the depth in DEPTHS
            idx = np.argmin(np.abs(DEPTHS - d))
            c = plt.get_cmap("tab20")(idx)

            if prevD is not None and d == prevD:
                idx_ls += 1
            else:
                idx_ls = 0
                prevD = d

            idx_tmp = (
                (idx_ls % len(LINES))
                if not plotBFSandFFSdata
                else ((idx_ls % len(LINES) - 1) + 1)
            )  # skip the first line style for the first case of each depth, as it is already used for the BFS case
            ls = LINES[
                idx_tmp
            ]  # skip the first line style for the first case of each depth, as it is already used for the BFS case
            addNxplot(d, w, x, Nx, ax, c, ls)
            plotd_vs_dNx(d, w, dN, ax2, plotBFS=False)

    depths = np.asarray(depths)
    widths = np.asarray(widths)
    reynolds = np.asarray(reynolds)
    deltaNs = np.asarray(deltaNs)

    plotd_vs_dNx(0, 0, 0, ax2, plotBFS=True)

    # GPR interpolation and plotting (separated)
    gpr_result = gpr_interpolateWithRe(
        depths, widths, reynolds, deltaNs, target_reynolds=re
    )
    write_deltaN_file(gpr_result, filenameDeltaN)
    plot_deltaN_grid(gpr_result)

    ax.set_title(f"n(x) factor using {'Loo' if doLoo else 'L2'} norm")
    ax.set_xlabel("x")
    ax.set_ylabel("n(x)")
    ax.grid()
    ax.legend(loc="center right", bbox_to_anchor=(1., 0.5), borderaxespad=0)
    ax2.set_title("Delta N as function of d for different w")
    ax2.set_xlabel("d")
    ax2.set_ylabel("Delta N")
    ax2.grid()
    ax2.legend()

    plt.show()


if __name__ == "__main__":
    main()
