import matplotlib.pyplot as plt
import numpy as np
from pp.colors import colors
from pp.fileManagement import extract_depth_width, readErrorHistoryMultiple
import os

plt.style.use("plots/style/jfm.mplstyle")


def main():
    """
    Main function to plot history points from one or multiple folders as a function of time or in phase space (u vs v). The points are stored in the historyPoints.dat file in each folder.
    """

    dir = "/home/victor/Desktop/PhD/data/periodDoubling/"
    script_path = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(script_path, "../../../images/periodDoubling.pdf")
    folders = ["d3_w33", "d3_w33.5"]
    folders = np.array([dir + f for f in folders])

    data = readErrorHistoryMultiple(folders)

    latex_width_cm = 8.0
    fig_width_in = latex_width_cm / 2.54

    _, ax = plt.subplots(
        figsize=(fig_width_in, fig_width_in * 0.75),
        constrained_layout=False,
    )

    col = ["tab:blue", "tab:orange"]
    line_styles = ["-", "--"]
    widths = [33, 33.5]
    periods = [280, 560]

    for i, f in enumerate(folders):
        _, errors = data[f]

        uL2, vL2 = errors[:, 0], errors[:, 3]
        vL2 = vL2[:periods[i]]
        uL2 = uL2[:periods[i]]

        # standarize the data by subtracting the mean and dividing by the standard deviation
        vL2 = (vL2 - np.mean(vL2)) / np.std(vL2)
        uL2 = (uL2 - np.mean(uL2)) / np.std(uL2)

        ax.plot(
            uL2,
            vL2,
            color=col[i],
            linestyle=line_styles[i],
            label=rf"$w/\delta^* = {widths[i]}$"
        )

    ax.set_xlabel(r"$\|u\|_2^\text{s}$")
    ax.set_ylabel(r"$\|v\|_2^\text{s}$", labelpad=10, rotation=0)
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-1, 0, 1])

    plt.legend()

    plt.savefig(save_path, format="pdf", bbox_inches="tight")
    print(colors.OKGREEN + f"✓ Plot saved to: {save_path}" + colors.ENDC)


if __name__ == "__main__":
    main()
