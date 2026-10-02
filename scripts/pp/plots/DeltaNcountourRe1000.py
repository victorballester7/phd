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
import matplotlib as mpl
from pp.fileManagement import extract_depth_width
from pp.hopfDataPoints import hopf2d_re_stable, hopf_re_unstable
from pp.filterData import average_curves
from plots.deltaNstyle import (
    COL_EXPERIMENT,
    add_excluded_region,
    add_numeric_contours,
    deltaN_colormap,
)
from pp.figureFrame import FigureFrame

plt.style.use("plots/style/jfm.mplstyle")

def addBifurcationMask(ax):
    """Hatch the region above the bifurcation curve, where no equilibrium exists."""
    add_excluded_region(ax, 1000, label="2D stability boundary")
    return


def experimental_levelsets_to_grid(w_range=(0, 90), d_range=(0, 4), resolution=200):
    """
    Convert experimental level set data to a grid where Z values are category integers.

    The function assigns each grid point to the category of the nearest level set,
    effectively creating regions bounded by the experimental contours.

    Parameters:
        experimentalDataFile: Path to CSV with columns [width, depth, category]
        w_range: (min, max) for width axis
        d_range: (min, max) for depth axis
        resolution: Number of grid points in each dimension

    Returns:
        X_grid, Y_grid, Z_grid: Meshgrid arrays where Z contains integer categories 1-4
    """
    # Create regular grid
    w_grid = np.linspace(w_range[0], w_range[1], resolution)
    d_grid = np.linspace(d_range[0], d_range[1], resolution)
    X_grid, Y_grid = np.meshgrid(w_grid, d_grid)

    # approximate model:
    # deltaN = 0.1 * w * tanh(44 *d/ w)
    Z_grid = 0.1 * X_grid * np.tanh(44 * Y_grid / X_grid)

    return X_grid, Y_grid, Z_grid


def addExperimentalContours(ax: plt.Axes, color=COL_EXPERIMENT):
    """
    Plot experimental data as filled contours matching the style of computed contours.
    """
    X_grid, Y_grid, Z_grid = experimental_levelsets_to_grid()

    levels = [1, 2, 3, 4, 5]  # Contour levels corresponding to categories

    # Draw contour lines at the boundaries
    CS = ax.contour(
        X_grid,
        Y_grid,
        Z_grid,
        levels=levels,
        colors=color,
        linestyles="--",
        linewidths=0.9,
        zorder=3,
    )

    ax.plot(
        [], [], color=color, linestyle="--", label="Experimental $\\Delta N$ levels\n(Crouch $\\it{et\\ al.}$ 2022)"
    )

    # Add labels at the center of each contour level
    # ax.clabel(CS, levels=levels, inline=True, fmt="%d")
    ax.clabel(
        CS,
        inline=True,
        manual=[(15, 0.5), (26, 0.5), (40, 0.75), (60, 1.2), (85, 1.5)],
        fmt="%d",
    )

    return CS


def plot_deltaN_report(xgrid, ygrid, zgrid, save_path: str):
    """Create a paper-quality figure for the Delta N contour plot."""
    # LaTeX document settings: 10pt font, figure width 6cm
    frame = FigureFrame(frame_w=7.0, aspect_ratio=0.8, pad_l=1.4, pad_b=0.8, pad_t=0.11)
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
        # cmap="Oranges_r",  # reversed: blue=low, white=high
        edgecolors="none",  # 🔴 remove cell edges
        linewidth=0,  # 🔴 ensure no lines
        antialiased=False,  # 🔴 avoids rendering seams
    )

    c.set_rasterized(True)

    # Contour lines (also masked)
    contour_levels = [1, 2, 3, 4, 5]
    add_numeric_contours(
        ax,
        xgrid,
        ygrid,
        zgrid,
        contour_levels,
        [(21, 0.5), (33, 0.7), (48, 0.95), (53, 1), (67, 1.2)],
    )
    ax.plot([], [], color="black", label=r"Numerical $\Delta n$ levels")  # for legend

    ax.set_xlim(0, 90)
    ax.set_ylim(0, 4)
    ax.set_yticks([0, 1, 2, 3, 4])

    # Add experimental data as contours (new method with proper formatting)
    addExperimentalContours(ax)
    addBifurcationMask(ax)

    # No grid: the shading already carries the value, and a grid over the
    # hatched region would be read as part of the hatch
    ax.grid(False)

    # draw ticks and spines above the colormap, contours and hatched region
    frame.axis_on_top()
    frame.jfm_ticks()

    # inside the hatched region: there is no data there, and putting the
    # legend outside made the panel half legend by width
    ax.legend(loc="upper right", frameon=False)

    # impose max value for colorbar is 6 to avoid too light colors for high values
    c.set_clim(0, 6)

    cbar = frame.colorbar(c)
    cbar.set_label(r"$\Delta n$", rotation=0, labelpad=10)
    ax.set_xlabel(r"$w$")
    ax.set_ylabel(r"$d$", rotation=0, labelpad=10)

    fig.savefig(save_path, format="pdf")
    print(colors.OKGREEN + f"Figure saved to {save_path}" + colors.ENDC)

    return fig


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


def main():
    """
    Computes Delta N for gap configurations and creates a report-quality contour plot
    using GPR interpolation.
    """
    n = 600
    doLoo = False
    re = 1000
    useOnly1Re = True 

    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))

    filenameDeltaN = f"../../../data/deltaN/DeltaN{re}.dat"
    filenameDeltaN = os.path.join(pathCurrentScript, filenameDeltaN)

    figurePath = "../../../images/DeltaNcountour.pdf"
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
    plot_deltaN_report(xgrid, ygrid, zgrid, figurePath)


if __name__ == "__main__":
    main()
