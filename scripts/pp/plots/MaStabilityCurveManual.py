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
    save_path = os.path.join(script_path, "../../../images/Mabifurcations.pdf")

    hopf2d = np.array(
        [
            [4, 17.5],
            [3.75, 18.5],
            [3.5, 20],
            [3.25, 20.5],
            [3, 25.4],
            [2.75, 26.5],
            [2.5, 28.5],
            [2.25, 31.5],
            [2, 39],
            [1.75, 46],
            [1.5, 59],
            [1.375, 65],
            [1.375, 70],
            [1.375, 75],
            [1.375, 90],
            [1.5, 95],
            [1.625, 105],
            [1.625, 120],
        ]
    )

    ma005 = np.array(
        [
            [4.0, 18.5],
            [3.5, 20.0],
            [3.0, 25.5],
            [2.5, 28.5],
            [2.0, 39.0],
            [1.5, 57.5],
        ]
    )

    ma02 = np.array(
        [
            [4.0, 16.5],
            [3.5, 17.5],
            [3.0, 18.5],
            [2.5, 24.5],
            [2.0, 29.0],
            [1.5, 37.5],
            [1.25, 50.0],
            [1.125, 60.0],
        ]
    )

    ma04 = np.array(
        [
            [4.0, 14.49],
            [3.5, 14.5],
            [3.0, 15.0],
            [2.5, 16.5],
            [2.0, 21.0],
            [1.5, 27.0],
            [1.0, 57.5],
            [0.875, 80.0],
        ]
    )

    ma06 = np.array(
        [
            [4.0, 12.5],
            [3.5, 13.49],
            [3.0, 13.5],
            [2.5, 14.5],
            [2.0, 19.0],
            [1.5, 27.5],
            [1.25, 35.5],
            [1.0, 52.5],
            [0.875, 60.0],
            [0.875, 80.0],
            [0.875, 120.0],
        ]
    )

    ma08 = np.array(
        [
            [4.0, 12.49],
            [3.5, 12.5],
            [3.0, 13.5],
            [2.5, 15.5],
            [2.0, 23.0],
            [1.5, 48.0],
            [1.25, 52.5],
            [1.125, 65.0],
            [1.125, 85.0],
        ]
    )

     # LaTeX document settings: 10pt font, figure width 6cm
    latex_width_cm = 7.85 # linewidth of tsfp template
    latex_font_pt = 9 # 9pt small font size in tsfp template (figures captions)
    
    fig_width_in, fig_height_in = compute_figure_size(latex_width_cm, aspect_ratio=0.75)

    # we make the figure as if there were only one column
    fig_width_in *= 2
    fig_height_in *= 2
    print(f"Figure size (inches): {fig_width_in:.2f} x {fig_height_in:.2f}")
    compute_mpl_fontsize(fig_width_in=fig_width_in, latex_width_cm=latex_width_cm, latex_font_pt=latex_font_pt)
    fig, ax = plt.subplots(figsize=(fig_width_in, fig_height_in), constrained_layout=True)

    x_hopf_smooth, y_hopf_smooth = smooth_curve(hopf2d[:, 1], hopf2d[:, 0])
    x_ma05_smooth, y_ma05_smooth = smooth_curve(ma005[:, 1], ma005[:, 0])
    x_ma2_smooth, y_ma2_smooth = smooth_curve(ma02[:, 1], ma02[:, 0])
    x_ma4_smooth, y_ma4_smooth = smooth_curve(ma04[:, 1], ma04[:, 0])
    x_ma6_smooth, y_ma6_smooth = smooth_curve(ma06[:, 1], ma06[:, 0])
    x_ma8_smooth, y_ma8_smooth = smooth_curve(ma08[:, 1], ma08[:, 0])

    # x_hopf_smooth, y_hopf_smooth = hopf2d[:, 1], hopf2d[:, 0]
    # x_ma05_smooth, y_ma05_smooth = ma005[:, 1], ma005[:, 0]
    # x_ma2_smooth, y_ma2_smooth = ma02[:, 1], ma02[:, 0]
    # x_ma4_smooth, y_ma4_smooth = ma04[:, 1], ma04[:, 0]
    # x_ma6_smooth, y_ma6_smooth = ma06[:, 1], ma06[:, 0]
    # x_ma8_smooth, y_ma8_smooth = ma08[:, 1], ma08[:, 0]

    ax.plot(
        x_hopf_smooth,
        y_hopf_smooth,
        linestyle="-",
        label="Inc. flow",
    )

    ax.plot(
        x_ma05_smooth,
        y_ma05_smooth,
        linestyle="-.",
        label=r"$\textrm{Ma}=0.05$",
    )

    ax.plot(
        x_ma2_smooth,
        y_ma2_smooth,
        linestyle="--",
        label=r"$\textrm{Ma}=0.2$",
    )

    ax.plot(
        x_ma4_smooth,
        y_ma4_smooth,
        linestyle=":",
        label=r"$\textrm{Ma}=0.4$",
    )

    ax.plot(
        x_ma6_smooth,
        y_ma6_smooth,
        linestyle=(0, (3, 1, 1, 1)),
        label=r"$\textrm{Ma}=0.6$",
    )

    ax.plot(
        x_ma8_smooth,
        y_ma8_smooth,
        linestyle=(0, (5, 5)),
        label=r"$\textrm{Ma}=0.8$",
    )

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

    # Tight layout
    plt.savefig(save_path, format="pdf")
    print(colors.OKGREEN + f"✓ Plot saved to: {save_path}" + colors.ENDC)


if __name__ == "__main__":
    main()
