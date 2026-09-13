import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from pp.DeltaN_computation import computeAmplitude, computeNx
from pp.fileManagement import extract_depth_width

RE_REFERENCE = 1000.0
BLASIUS_C = 1.7207876573

plt.style.use("plots/style/jfm.mplstyle")

linestyles = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1))]

n = 5
cmap_tmp = plt.get_cmap("Purples_r")
colors = cmap_tmp(np.linspace(0, 0.7, n))


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(script_path, "../../../images/L2POvsX.pdf")
    cases = {
        "d1.5_w60": "w=60",
        "d1.5_w63": "w=63",
        "d1.5_w65": "w=65",
        "d1.5_w70": "w=70",
        "d1.5_w80": "w=80",
    }

    latex_width_cm = 6.5
    fig_width_in = latex_width_cm / 2.54
    _, ax = plt.subplots(
        figsize=(fig_width_in, fig_width_in * 0.75),
        constrained_layout=False,
    )
    re = 1000
    ylabel = r"$\||\boldsymbol{u}'|\|_{L^2(\mathcal{Y}(x))}$"
    # ylabel = r"$\omega_r$"

    for i, case in enumerate(cases.keys()):
        file_path = f"/home/victor/Desktop/PhD/src/incGapRe{re}/baseflow/dns/{case}/data/pointsPO_n600.dat"
        x, L2 = computeAmplitude(file_path, False, field="|u|")
        # window average L2 norm to smooth out the curve
        L2_smooth = gaussian_filter1d(L2, sigma=10)
        # L2_smooth = np.convolve(L2, np.ones(window_size) / window_size, mode="same")

        col = colors[i]
        # lin = linestyles[i]
        lin = "-"
        _, width = extract_depth_width(case)
        x = x - width
        ax.plot(x, L2, linestyle="-", color=col, linewidth=0.75, alpha=0.3)
        ax.plot(
            x, L2_smooth, linestyle=lin, color=col, label=rf"${cases[case]}\delta^*$"
        )

    ax.set_xlabel(r"$(x - w)/\delta^*$")
    ax.set_xlim(0, 900)
    ax.set_ylabel(ylabel, rotation=0, labelpad=30)

    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
    )
    # plt.show()
    plt.savefig(
        save_path,
        format="pdf",
        bbox_inches="tight",
    )


if __name__ == "__main__":
    main()
