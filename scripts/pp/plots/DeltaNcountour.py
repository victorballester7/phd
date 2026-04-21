import numpy as np
import matplotlib.pyplot as plt
import os
from pp.codenamesNfactor import code_names
from pp.DeltaN_computation import (
    computeNx,
    computeDeltaN,
    gpr_interpolate,
    write_deltaN_file,
)
from pp.colors import colors
import matplotlib as mpl
import matplotlib.colors as mcolors
from pp.fontSizeLaTex import compute_mpl_fontsize, compute_figure_size
from pp.fileManagement import extract_width_depth

plt.style.use("plots/style/tsfp.mplstyle")

# Bifurcation curve data
BIFURCATION_CURVE = np.array(
    # [
    #     [16.25, 4.0],
    #     [15.75, 3.75],
    #     [15.25, 3.5],
    #     [15.25, 3.25],
    #     [15.25, 3.0],
    #     [15.25, 2.75],
    #     [20.5, 2.5],
    #     [31.5, 2.25],
    #     [39.0, 2.0],
    #     [46.0, 1.75],
    #     [56.5, 1.5],
    #     [59.0, 1.46],
    #     [65.0, 1.375],
    #     [70.0, 1.375],
    #     [75.0, 1.375],
    #     [80.0, 1.375],
    #     [100.0, 1.375],
    #     [120.0, 1.375],
    #     [130.0, 1.375],
    # ]
    [
        [17.5, 4],
        [18.5, 3.75],
        [20, 3.5],
        [20.5, 3.25],
        [25.4, 3],
        [26.5, 2.75],
        [28.5, 2.5],
        [31.5, 2.25],
        [39, 2],
        [46, 1.75],
        [59, 1.5],
        [65, 1.375],
        [70, 1.375],
        [75, 1.375],
        [80, 1.375],
        [100, 1.375],
        [120, 1.375],
        [130, 1.375],
    ]
)


def addBifurcationMask(ax):
    """Create a mask for points above/right of the bifurcation curve."""
    # Interpolate bifurcation curve: w as function of d
    x_curve = BIFURCATION_CURVE[:, 0]
    y_curve = BIFURCATION_CURVE[:, 1]

    # fill everything ABOVE the curve with white
    ax.fill_between(x_curve, y_curve, y2=plt.ylim()[1], color="white", zorder=2)

    ax.plot(
        x_curve,
        y_curve,
        color="black",
        linestyle=":",
        label="2D instability boundary",
        zorder=2,
    )
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


def addExperimentalContours(ax: plt.Axes, color="tab:Orange"):
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
        zorder=1,
    )

    ax.plot([], [], color=color, linestyle="--", label=r"Experimental $\Delta n$ levels")

    # Add labels at the center of each contour level
    # ax.clabel(CS, levels=levels, inline=True, fmt="%d")
    ax.clabel(
        CS, inline=True, manual=[(15, 0.5), (26, 0.5), (40, 0.75), (60, 1.2), (85, 1.5)], fmt="%d"
    )

    return CS


def plot_deltaN_report(result, save_path: str):
    """Create a paper-quality figure for the Delta N contour plot."""
    # LaTeX document settings: 10pt font, figure width 6cm
    latex_width_cm = 7.85  # linewidth of tsfp template
    latex_font_pt = 9  # 9pt small font size in tsfp template (figures captions)

    fig_width_in, fig_height_in = compute_figure_size(latex_width_cm, aspect_ratio=0.9)

    # we make the figure as if there were only one column
    fig_width_in *= 2
    fig_height_in *= 2
    print(f"Figure size (inches): {fig_width_in:.2f} x {fig_height_in:.2f}")
    compute_mpl_fontsize(
        fig_width_in=fig_width_in,
        latex_width_cm=latex_width_cm,
        latex_font_pt=latex_font_pt,
    )
    fig, ax = plt.subplots(
        figsize=(fig_width_in, fig_height_in), constrained_layout=True
    )

    # Apply bifurcation mask to Z_grid
    # mask = get_bifurcation_mask(result.X_grid, result.Y_grid)
    # Z_masked = np.ma.array(result.Z_grid, mask=mask)

    def truncate_colormap(cmap, minval=0.3, maxval=1.0, n=256):
        return mcolors.LinearSegmentedColormap.from_list(
            f"trunc({cmap.name},{minval:.2f},{maxval:.2f})",
            cmap(np.linspace(minval, maxval, n)),
        )

    cmap = truncate_colormap(plt.get_cmap("Blues_r"), 0.15, 1.0)

    c = ax.pcolormesh(
        result.X_grid,
        result.Y_grid,
        result.Z_grid,
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
    CS = ax.contour(
        result.X_grid,
        result.Y_grid,
        result.Z_grid,
        levels=contour_levels,
        colors="darkblue",
        zorder=2,
    )
    ax.plot([], [], color="darkblue", label=r"Computed $\Delta n$ levels")  # for legend
    ax.clabel(
        CS,
        inline=True,
        manual=[(21, 0.5), (33, 0.7), (48, 0.95), (53, 1), (67, 1.2)],
        fmt="%d",
    )

    # Add experimental data as contours (new method with proper formatting)
    addExperimentalContours(ax)
    addBifurcationMask(ax)

    # Draw grid lines manually on top of everything (including white mask)
    ax.set_axisbelow(False)

    # ax.legend(loc="upper right")
    # add legend below the plot, centered
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2))

    # impose max value for colorbar is 6 to avoid too light colors for high values
    c.set_clim(0, 6)

    cbar = fig.colorbar(c, ax=ax, orientation="vertical",aspect=15)
    cbar.set_label(r"$\Delta n$", rotation=0, labelpad=15)
    ax.set_xlabel(r"$w/\delta^*$")
    ax.set_ylabel(r"$d/\delta^*$", rotation=0, labelpad=25)
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

    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))

    filenameDeltaN = "../../../latex/Images/data/DeltaN.dat"
    filenameDeltaN = os.path.join(pathCurrentScript, filenameDeltaN)

    figurePath = "../../../images/DeltaNcountour.pdf"
    figurePath = os.path.join(pathCurrentScript, figurePath)

    basePath = "../../../src/incGapRe1000/directLinearSolver/blowingSuction/"
    basePath = os.path.join(pathCurrentScript, basePath)

    # Flat plate reference
    dataFile_flat = f"../../../src/flatPlateRe1000inc/directLinearSolver/blowingSuction/wgnInsideDomainDivFree/data/pointsavg_n{n}.dat"
    dataFile_flat = os.path.join(pathCurrentScript, dataFile_flat)
    x_flat, Nx_flat = computeNx(dataFile_flat, doLoo)

    depths = np.array([])
    widths = np.array([])
    deltaNs = np.array([])

    for dw in code_names:
        dataFile_dw = os.path.join(basePath, dw, "data", f"pointsavg_n{n}.dat")
        d, w = extract_width_depth(dataFile_dw)
        x, Nx = computeNx(dataFile_dw, doLoo)
        depths = np.append(depths, d)
        widths = np.append(widths, w)
        dN = computeDeltaN(w, x, Nx, x_flat, Nx_flat)
        deltaNs = np.append(deltaNs, dN)

    # GPR interpolation
    gpr_result = gpr_interpolate(depths, widths, deltaNs,w_max=90)
    write_deltaN_file(gpr_result, filenameDeltaN)

    # Create report figure
    plot_deltaN_report(gpr_result, figurePath)


if __name__ == "__main__":
    main()
