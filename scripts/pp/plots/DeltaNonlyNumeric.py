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
import matplotlib.colors as mcolors
from pp.fileManagement import extract_depth_width
from pp.hopfDataPoints import hopf2d_re_stable, hopf_re_unstable
from pp.filterData import average_curves

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

    def __init__(self, re):
        self.re = re
        self.colorCountour = ""
        self.colorLevelSet = ""
        self.contourLevels = []
        self.labelPositions = []
        match re:
            case 800:
                self.colorCountour = "Greens_r"
                self.colorLevelSet = "darkgreen"
                self.contourLevels = [1, 2, 3, 4, 5, 6]
                self.labelPositions = [(21, 0.5), (33, 0.7), (44, 0.95), (53, 1), (67, 1.2), (82, 1.4)]
            case 1000:
                self.colorCountour = "Blues_r"
                self.colorLevelSet = "darkblue"
                self.contourLevels = [1, 2, 3, 4, 5]
                self.labelPositions = [(21, 0.5), (33, 0.7), (44, 0.95), (53, 1), (67, 1.2)]
            case 3000:
                self.colorCountour = "Reds_r"
                self.colorLevelSet = "maroon"
                self.contourLevels = [1, 2, 3, 4]
                self.labelPositions = [(21, 0.5), (33, 0.7), (48, 0.8), (63, 0.75)]

def addBifurcationMask(ax, re):
    """Create a mask for points above/right of the bifurcation curve."""
    # Interpolate bifurcation curve: w as function of d
    
    bif_curve = average_curves(hopf2d_re_stable[re], hopf_re_unstable[re], npts=500)[2]

    x_curve = bif_curve[:, 0]
    y_curve = bif_curve[:, 1]

    # fill everything ABOVE the curve with white
    ax.fill_between(x_curve, y_curve, y2=plt.ylim()[1], color="white", zorder=2)

    ax.plot(
        x_curve,
        y_curve,
        color="black",
        linestyle=":",
        label="Linear stability boundary",
        zorder=2,
    )
    return


def plot_deltaN_report(restyle, xgrid, ygrid, zgrid, save_path: str):
    """Create a paper-quality figure for the Delta N contour plot."""
    # LaTeX document settings: 10pt font, figure width 6cm
    latex_width_cm = 6.5
    fig_width_in = latex_width_cm / 2.54

    fig, ax = plt.subplots(
        figsize=(fig_width_in, fig_width_in * 0.75),
        constrained_layout=False,
    )

    # Apply bifurcation mask to Z_grid
    # mask = get_bifurcation_mask(result.X_grid, result.Y_grid)
    # Z_masked = np.ma.array(result.Z_grid, mask=mask)

    def truncate_colormap(cmap, minval=0.3, maxval=1.0, n=256):
        return mcolors.LinearSegmentedColormap.from_list(
            f"trunc({cmap.name},{minval:.2f},{maxval:.2f})",
            cmap(np.linspace(minval, maxval, n)),
        )


    cmap = truncate_colormap(plt.get_cmap(restyle.colorCountour), 0.15, 1.0)

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


    CS = ax.contour(
        xgrid,
        ygrid,
        zgrid,
        levels=restyle.contourLevels,
        colors=restyle.colorLevelSet,
        zorder=2,
    )

    ax.clabel(
        CS,
        inline=True,
        manual=restyle.labelPositions,
        fmt="%d",
    )

    # Add experimental data as contours (new method with proper formatting)
    addBifurcationMask(ax, restyle.re)

    # Draw grid lines manually on top of everything (including white mask)
    ax.set_axisbelow(False)

    # ax.legend(loc="upper right")
    # add legend below the plot, centered
    # ax.legend(loc="center left", bbox_to_anchor=(1.3, 0.5))
    # ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3))

    # impose max value for colorbar is 6 to avoid too light colors for high values
    c.set_clim(0, 7 if restyle.re == 800 else 5)

    cbar = fig.colorbar(c, ax=ax, orientation="vertical", aspect=15)
    cbar.set_ticks([0, 1, 2, 3, 4, 5, 6, 7] if restyle.re == 800 else [0, 1, 2, 3, 4, 5])

    # set min max values for colorbar colormap

    cbar.set_label(r"$\Delta n$", rotation=0, labelpad=10)
    ax.set_xlabel(r"$w/\delta^*$")
    ax.set_ylabel(r"$d/\delta^*$", rotation=0, labelpad=10)
    ax.set_xticks([0, 20, 40, 60, 80])
    ax.set_xlim(0, 90)
    ax.set_ylim(0, 4)

    # Set white background for masked region

    fig.savefig(save_path, format="pdf")
    print(colors.OKGREEN + f"Figure saved to {save_path}" + colors.ENDC)

    return fig


def main():
    """
    Computes Delta N for gap configurations and creates a report-quality contour plot
    using GPR interpolation.
    """
    n = 600
    doLoo = False
    re = 800
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
