import numpy as np
import matplotlib.pyplot as plt
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pp.codenamesNfactor import code_names_Re
from pp.DeltaN_computation import (
    computeNx,
    computeDeltaN,
    gpr_interpolateWithRe,
    write_deltaN_file,
)
from pp.colors import colors
from pp.fileManagement import extract_depth_width
from pp.hopfDataPoints import hopf2d_re_stable, hopf_re_unstable
from pp.filterData import average_curves
from plots.deltaNstyle import add_excluded_region, add_numeric_contours, deltaN_colormap
from pp.figureFrame import FigureFrame

plt.style.use("plots/style/jfm.mplstyle")


_worker_x_flat = None
_worker_Nx_flat = None


def init_delta_n_worker(x_flat, Nx_flat):
    global _worker_x_flat, _worker_Nx_flat
    _worker_x_flat = x_flat
    _worker_Nx_flat = Nx_flat


def compute_delta_n_case(data_file, doLoo):
    """Load one gap case and compute its Delta N independently."""
    d, w = extract_depth_width(data_file)
    x, Nx = computeNx(data_file, doLoo)
    dN = computeDeltaN(w, x, Nx, _worker_x_flat, _worker_Nx_flat)
    return d, w, dN


class ReStyle:

    # The colour ramp is shared with figure DeltaNcountour.pdf (see
    # plots.deltaNstyle): three Reynolds numbers drawn in three different hues
    # cannot be compared by eye, which is the whole point of putting them side
    # by side. Only the contour levels and their label positions differ.
    def __init__(self, re):
        self.re = re
        self.contourLevels = []
        self.labelPositions = []
        match re:
            case 800:
                self.contourLevels = [1, 2, 3, 4, 5, 6]
                self.labelPositions = [(21, 0.5), (33, 0.7), (44, 0.95), (53, 1), (67, 1.2), (82, 1.4)]
            case 1000:
                self.contourLevels = [1, 2, 3, 4, 5]
                self.labelPositions = [(21, 0.5), (33, 0.7), (44, 0.95), (53, 1), (67, 1.2)]
            case 3000:
                self.contourLevels = [1, 2, 3, 4]
                self.labelPositions = [(21, 0.5), (33, 0.7), (48, 0.8), (63, 0.75)]

def addBifurcationMask(ax, re):
    """Hatch the region above the bifurcation curve, where no equilibrium exists."""
    add_excluded_region(ax, re, label="Linear stability boundary")
    return


def plot_deltaN_report(restyle, xgrid, ygrid, zgrid, save_path: str):
    """Create a paper-quality figure for the Delta N contour plot."""
    # LaTeX document settings: 10pt font, figure width 6cm
    frame = FigureFrame(frame_w=4.1, aspect_ratio=0.8, pad_l=1.3, pad_b=0.8, pad_t=0.11)
    fig, ax = frame.fig, frame.ax

    # Apply bifurcation mask to Z_grid
    # mask = get_bifurcation_mask(result.X_grid, result.Y_grid)
    # Z_masked = np.ma.array(result.Z_grid, mask=mask)

    cmap = deltaN_colormap()

    c = ax.pcolormesh(
        xgrid,
        ygrid,
        zgrid,
        shading="auto",
        cmap=cmap,
        edgecolors="none",  # 🔴 remove cell edges
        linewidth=0,  # 🔴 ensure no lines
        antialiased=False,  # 🔴 avoids rendering seams
    )

    c.set_rasterized(True)


    ax.set_xlim(0, 90)
    ax.set_ylim(0, 4)

    add_numeric_contours(
        ax,
        xgrid,
        ygrid,
        zgrid,
        restyle.contourLevels,
        restyle.labelPositions,
    )

    addBifurcationMask(ax, restyle.re)

    # No grid: the shading already carries the value, and a grid over the
    # hatched region would be read as part of the hatch
    ax.grid(False)

    # draw ticks and spines above the colormap, contours and hatched region
    frame.axis_on_top()
    frame.jfm_ticks()

    # ax.legend(loc="upper right")
    # add legend below the plot, centered
    # ax.legend(loc="center left", bbox_to_anchor=(1.3, 0.5))
    # ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3))

    # impose max value for colorbar is 6 to avoid too light colors for high values
    c.set_clim(0, 7)
    # c.set_clim(0, 7 if restyle.re == 800 else 5)

    cbar = frame.colorbar(c)
    # cbar.set_ticks([0, 1, 2, 3, 4, 5, 6, 7] if restyle.re == 800 else [0, 1, 2, 3, 4, 5])
    cbar.set_ticks([0, 1, 2, 3, 4, 5, 6, 7])

    # set min max values for colorbar colormap

    cbar.set_label(r"$\Delta n$", rotation=0, labelpad=10)
    ax.set_xlabel(r"$w$")
    ax.set_ylabel(r"$d$", rotation=0, labelpad=10)
    ax.set_xticks([0, 20, 40, 60, 80])

    fig.savefig(save_path, format="pdf")
    print(colors.OKGREEN + f"Figure saved to {save_path}" + colors.ENDC)

    return fig


def main():
    """Both panels of the Delta n figure, at Re = 800 and Re = 3000."""
    for re in (800, 3000):
        plot_one_reynolds(re)


def plot_one_reynolds(re):
    """
    Computes Delta N for gap configurations and creates a report-quality contour plot
    using GPR interpolation.
    """
    n = 600
    doLoo = False
    useOnly1Re = True
    restyle = ReStyle(re)

    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))

    filenameDeltaN = f"../../../data/deltaN/DeltaN{re}.dat"
    filenameDeltaN = os.path.join(pathCurrentScript, filenameDeltaN)

    figurePath = f"../../../images/DeltaNcountourOnlyNumericRe{re}.pdf"
    figurePath = os.path.join(pathCurrentScript, figurePath)

    basePath = "../../../src/incGapRe1000/directLinearSolver/blowingSuction/"
    basePath = os.path.join(pathCurrentScript, basePath)

    # Flat plate reference
    dataFile_flat = f"../../../src/flatPlateRe1000inc/directLinearSolver/blowingSuction/data/pointsavg_n{n}.dat"
    dataFile_flat = os.path.join(pathCurrentScript, dataFile_flat)
    x_flat, Nx_flat = computeNx(dataFile_flat, doLoo)

    depths = []
    widths = []
    reynolds = []
    deltaNs = []

    # check if DeltaN file already exists, if so, skip computation and just plot
    recompute = True
    if os.path.exists(filenameDeltaN):
        print(
            colors.WARNING
            + f"DeltaN file {filenameDeltaN} already exists. Would you like to use the current file (y) or recompute the whole simulations (n)? (Y/n): "
            + colors.ENDC,
            end="",
        )
        user_input = input().strip().lower()
        if user_input == "y" or user_input == "":
            print(
                colors.OKGREEN
                + f"Using existing DeltaN file {filenameDeltaN}"
                + colors.ENDC
            )

            recompute = False

    if recompute:
        for r in code_names_Re.keys():
            if useOnly1Re and r != re:
                continue
            flat_plate_file = (
                f"../../../src/flatPlateRe{r}inc/"
                "directLinearSolver/blowingSuction/"
                f"data/pointsavg_n{n}.dat"
            )

            flat_plate_file = os.path.join(
                pathCurrentScript,
                flat_plate_file,
            )

            x_flat, Nx_flat = computeNx(flat_plate_file, doLoo)

            basePath = f"../../../src/incGapRe{r}/directLinearSolver/blowingSuction/"

            basePath = os.path.join(
                pathCurrentScript,
                basePath,
            )

            data_files = [
                os.path.join(basePath, dw, "data", f"pointsavg_n{n}.dat")
                for dw in code_names_Re[r]
            ]
            with ProcessPoolExecutor(
                max_workers=min(len(data_files), os.cpu_count() or 1),
                initializer=init_delta_n_worker,
                initargs=(x_flat, Nx_flat),
            ) as executor:
                results = executor.map(
                    compute_delta_n_case,
                    data_files,
                    repeat(doLoo),
                    chunksize=1,
                )

            for d, w, dN in results:
                depths.append(d)
                widths.append(w)
                reynolds.append(r)
                deltaNs.append(dN)

        depths = np.asarray(depths)
        widths = np.asarray(widths)
        reynolds = np.asarray(reynolds)
        deltaNs = np.asarray(deltaNs)

        # GPR interpolation
        gpr_result = gpr_interpolateWithRe(
            depths,
            widths,
            reynolds,
            deltaNs,
            target_reynolds=re,
        )
        write_deltaN_file(gpr_result, filenameDeltaN)

        xgrid, ygrid, zgrid = gpr_result.X_grid, gpr_result.Y_grid, gpr_result.Z_grid
    else:
        # Load data from file
        data = np.loadtxt(filenameDeltaN)
        x = data[:, 0]
        y = data[:, 1]
        z = data[:, 2]
        

        n = int(np.sqrt(len(z)))

        xgrid = x.reshape(n, n)
        ygrid = y.reshape(n, n)
        zgrid = z.reshape(n, n)
    # Create report figure
    plot_deltaN_report(restyle, xgrid, ygrid, zgrid, figurePath)


if __name__ == "__main__":
    main()
