import numpy as np
import matplotlib.pyplot as plt
import os
from pp.colors import mix_with_black, colors
from scipy.interpolate import PchipInterpolator
import matplotlib as mpl
from scipy.interpolate import make_interp_spline
from scipy.interpolate import interp1d
from pp.fontSizeLaTex import compute_mpl_fontsize, compute_figure_size


# plt.style.use("plots/style/tsfp.mplstyle")
plt.style.use('plots/style/customvictor.mplstyle')


def smooth_curve(x, y, num_points=200):
    # Sort x and y based on x to ensure strictly increasing for PCHIP
    sort_idx = np.argsort(x)
    x_sorted = x[sort_idx]
    y_sorted = y[sort_idx]

    # Create PCHIP interpolator
    pchip = PchipInterpolator(x_sorted, y_sorted)

    # Generate new x values for smooth curve
    x_smooth = np.linspace(x_sorted.min(), x_sorted.max(), num_points)
    y_smooth = pchip(x_smooth)

    return x_smooth, y_smooth


def extend_2d_curve(x_2d, y_2d, x_3d, y_3d):
    idx = 0
    for i in range(len(x_3d)):
        x = x_3d[i]
        if x < 40:
            continue

        # interpolate y value at this x using the 2D curve
        if x < x_2d.min() or x > x_2d.max():
            continue
        y_interp = np.interp(x, x_2d, y_2d)

        if y_interp < y_3d[i]:
            idx = i
            break

    x_3d = x_3d[:idx]
    y_3d = y_3d[:idx]

    # add the tail of the 2D curve to the 3D curve
    x_3d_last = x_3d[-1]
    x_3d = np.concatenate((x_3d, x_2d[x_2d > x_3d_last]))
    y_3d = np.concatenate((y_3d, y_2d[x_2d > x_3d_last]))

    return x_3d, y_3d


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(script_path, "../../../images/3d2dbifurcations.pdf")

    hopf2d = np.array(
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
        ]
    )

    inst3d = np.array(
        [
            [7.5, 4],
            [7.8, 3.5],
            [8.3, 3],
            [8.7, 2.5],
            [9.4, 2.25],
            [22, 2],
            [30, 1.75],
            [54, 1.675],
        ]
    )

    # LaTeX document settings: 10pt font, figure width 6cm
    latex_width_cm = 7.85  # linewidth of tsfp template
    latex_font_pt = 9  # 9pt small font size in tsfp template (figures captions)

    fig_width_in, fig_height_in = compute_figure_size(latex_width_cm, aspect_ratio=0.75)

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

    x_hopf_smooth, y_hopf_smooth = smooth_curve(hopf2d[:, 0], hopf2d[:, 1])
    x_inst_smooth, y_inst_smooth = smooth_curve(inst3d[:, 0], inst3d[:, 1])

    x_inst_smooth_extended, y_inst_smooth_extended = extend_2d_curve(
        x_hopf_smooth, y_hopf_smooth, x_inst_smooth, y_inst_smooth
    )

    print(x_inst_smooth_extended.shape, y_inst_smooth_extended.shape)
    print(x_hopf_smooth.shape, y_hopf_smooth.shape)

    ax.plot(
        x_hopf_smooth,
        y_hopf_smooth,
        linestyle="--",
        label="2D Hopf bifurcation",
        zorder=3,
    )

    ax.plot(
        x_inst_smooth_extended,
        y_inst_smooth_extended,
        linestyle="-",
        label="3D instability boundary",
        zorder=2,
    )

    # ax.plot([19], [4], "* ", color="tab:purple", markersize=15)

    # smooth_inst_x, smooth_inst_y = smooth_curve(inst3d[:, 1], inst3d[:, 0])
    # ax.plot(smooth_inst_x, smooth_inst_y, linestyle='-', label='3D Instability', zorder=3)
    # ax.plot(
    #     inst3d[:, 1],
    #     inst3d[:, 0],
    #     "x-",
    #     color="blue",
    #     label="3D Instability Data Points",
    #     zorder=4,
    # )

    ax.set_xlabel(r"$w/\delta^*$")
    ax.set_ylabel(r"$d/\delta^*$", rotation=0, labelpad=20)

    # Grid and legend
    ax.legend()

    # Save with fixed figure size (no bbox_inches="tight" to preserve exact dimensions)
    plt.savefig(save_path, format="pdf")
    print(colors.OKGREEN + f"✓ Plot saved to: {save_path}" + colors.ENDC)


if __name__ == "__main__":
    main()
