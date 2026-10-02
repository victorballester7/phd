import matplotlib.pyplot as plt
import numpy as np
from pp.colors import colors as c
from pp.fileManagement import extract_depth_width, readErrorHistoryMultiple
import os
from pp.figureFrame import FigureFrame

plt.style.use("plots/style/jfm.mplstyle")


def add_direction_arrow(ax, x, y, frac, color):
    """A single arrowhead on the orbit, so the sense of travel is readable."""
    i = int(frac * (len(x) - 2))
    ax.annotate(
        "",
        xy=(x[i + 1], y[i + 1]),
        xytext=(x[i], y[i]),
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            linewidth=0.0,
            mutation_scale=9,
            shrinkA=0,
            shrinkB=0,
        ),
    )


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

    frame = FigureFrame(frame_w=4.5, aspect_ratio=1.0, pad_l=2.4, pad_b=0.9, pad_t=0.05)
    fig, ax = frame.fig, frame.ax

    cmap = plt.get_cmap('Dark2')
    col = [cmap.colors[0], cmap.colors[6]]  
    line_styles = ["-", "--"]
    widths = [33.0, 33.5]
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
            label=rf"$w = {widths[i]}$"
        )
        for frac in (0.15, 0.65):
            add_direction_arrow(ax, uL2, vL2, frac, col[i])

    ax.set_xlabel(r"$\|u\|_2^\text{s}$")
    ax.set_ylabel(r"$\|v\|_2^\text{s}$", labelpad=10, rotation=0)
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-1, 0, 1])
    # equal spans on both axes, with a margin around the orbit
    ax.set_xlim(-2.05, 1.5)
    ax.set_ylim(-1.5, 2.05)
    # both axes are the same standardised norm, so the orbit should not be
    # stretched in one direction; "datalim" keeps the square frame fixed and
    # adjusts the limits instead of shrinking the axes box
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(False)
    frame.axis_on_top()
    frame.jfm_ticks()

    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)

    fig.savefig(save_path, format="pdf")
    print(c.OKGREEN + f"✓ Plot saved to: {save_path}" + c.ENDC)


if __name__ == "__main__":
    main()
