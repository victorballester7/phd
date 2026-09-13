import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from matplotlib.ticker import FixedLocator, NullLocator

RE_REFERENCE = 1000.0
BLASIUS_C = 1.7207876573

plt.style.use("plots/style/jfm.mplstyle")

linestyles = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1))]

# cmap = plt.get_cmap("tab10")
col_fp = "tab:blue"
col_bfs = "tab:orange"



def read_neutral_curve(file_path):
    with open(file_path) as f:
        text = f.read().strip()

    # Split on one or more blank lines
    blocks = text.split("\n\n")

    arrays = []
    reynolds = []
    omegas = []
    x = []
    for i, block in enumerate(blocks):
        if not block.strip():
            raise ValueError(f"Block {i} is empty in file {file_path}")

        arr = np.loadtxt(block.splitlines())
        arrays.append(arr)

        reynolds.append(arr[:, 0])
        omegas.append(arr[:, 1])
        x.append(arr[:, 2])

    return reynolds, omegas, x


def main():
    depthConstant = True
    useReynolds = True

    if depthConstant:
        n = 6
        cmap_tmp = plt.get_cmap("Greens")
        colors = cmap_tmp(np.linspace(0.2, 1.0, n))
    else:
        n = 7
        cmap_tmp = plt.get_cmap("Reds")
        colors = cmap_tmp(np.linspace(0.2, 1.0, n))

    script_path = os.path.dirname(os.path.abspath(__file__))
    if depthConstant:
        cases = {
            "d1.5_w10": "w=10",
            "d1.5_w15": "w=15",
            "d1.5_w20": "w=20",
            # "d1.5_w25": "w=25",
            "d1.5_w30": "w=30",
            # "d1.5_w35": "w=35",
            "d1.5_w40": "w=40",
            # "d1.5_w45": "w=45",
            "d1.5_w50": "w=50",
        }
        save_path = os.path.join(script_path, "../../../images/neutrald1.5.pdf")
    else:
        cases = {
            "d0.25_w40": "d=0.25",
            "d0.5_w40": "d=0.5",
            "d0.75_w40": "d=0.75",
            "d1_w39": "d=1",
            "d1.25_w40": "d=1.25",
            "d1.5_w40": "d=1.5",
            "d1.75_w40": "d=1.75",
        }
        save_path = os.path.join(script_path, "../../../images/neutralw40.pdf")

    latex_width_cm = 6.5
    fig_width_in = latex_width_cm / 2.54
    _, ax = plt.subplots(
        figsize=(fig_width_in, fig_width_in * 0.75),
        constrained_layout=False,
    )
    re = 1000
    ylabel = r"$F$"
    # ylabel = r"$\displaystyle\frac{\omega}{u_\infty/\delta^*(x)}$"

    # flat plate case
    file_path_fp = "/home/victor/Desktop/PhD/data/neutralCurves/flat/neutral_curve_temporal_gaster.dat"
    reynolds_fp, omegas_fp, x_fp = read_neutral_curve(file_path_fp)

    lin_fp = linestyles[1]
    for i, (r, o, x) in enumerate(zip(reynolds_fp, omegas_fp, x_fp)):
        x_plot = r.copy() if useReynolds else x.copy()
        o_plot = o.copy() * 1e6 / (r.copy() * BLASIUS_C)
        x_plot *= BLASIUS_C if useReynolds else 1.0
        # o_plot = o
        # o_plot = o / np.sqrt(
        #     1 + x * BLASIUS_C**2 / RE_REFERENCE
        # )  # convert to local frequency
        if i == 0:
            ax.plot(x_plot, o_plot, linestyle=lin_fp, color=col_fp, label="Flat plate")
        else:
            ax.plot(x_plot, o_plot, linestyle=lin_fp, color=col_fp)

    if depthConstant:
        file_path_bfs = "/home/victor/Desktop/PhD/data/neutralCurves/bfs/neutral_curve_temporal_gaster.dat"
        reynolds_bfs, omegas_bfs, x_bfs = read_neutral_curve(file_path_bfs)
        lin_bfs = linestyles[2]
        for i, (r, o, x) in enumerate(zip(reynolds_bfs, omegas_bfs, x_bfs)):
            x_plot = r.copy() if useReynolds else x.copy()
            o_plot = o.copy() * 1e6 / (r.copy() * BLASIUS_C)
            x_plot *= BLASIUS_C if useReynolds else 1.0
            if i == 0:
                ax.plot(x_plot, o_plot, linestyle=lin_bfs, color=col_bfs, label="BFS")
            else:
                ax.plot(x_plot, o_plot, linestyle=lin_bfs, color=col_bfs)

    for i, case in enumerate(cases.keys()):
        file_path = f"/home/victor/Desktop/PhD/data/neutralCurves/{case}/neutral_curve_temporal_gaster.dat"
        col = colors[i]
        lin = linestyles[(i + 1) % len(linestyles)]
        lin = "-"

        reynolds, omegas, xs = read_neutral_curve(file_path)

        for j, (r, o, x) in enumerate(zip(reynolds, omegas, xs)):
            x_plot = r if useReynolds else x
            o_plot = o * 1e6 / (r * BLASIUS_C)
            x_plot *= BLASIUS_C if useReynolds else 1.0
            if j == 0:
                ax.plot(
                    x_plot,
                    o_plot,
                    linestyle=lin,
                    color=col,
                    label=rf"${cases[case]}\delta^*$",
                )
            else:
                ax.plot(x_plot, o_plot, linestyle=lin, color=col)

    ax.set_xlabel(r"$\mbox{\textit{Re}}$" if useReynolds else r"$x/\delta^*$")
    ax.set_ylabel(ylabel, rotation=0, labelpad=5)

    def re_to_x(re):
        return (
            ((np.asarray(re) / RE_REFERENCE) ** 2 - 1.0) * RE_REFERENCE / BLASIUS_C**2
        )

    def x_to_re(x):
        return RE_REFERENCE * np.sqrt(1.0 + BLASIUS_C**2 * np.asarray(x) / RE_REFERENCE)

    #
    ax2 = ax.secondary_xaxis(
        "top",
        functions=(re_to_x, x_to_re),
    )
    ax2.set_xlabel(r"$x/\delta^*$")

    # ------------------------------------------------------------------
    # Zoomed inset
    # ------------------------------------------------------------------
    axins = inset_axes(
        ax,
        width="50%",  # width = 30% of parent_bbox
        height="50%",  # height : 30%
        loc="upper right",
        # bbox_to_anchor=(100,50,120,120),  # position relative to parent
    )

    for r, o, x in zip(reynolds_fp, omegas_fp, x_fp):
        x_plot = r.copy() if useReynolds else x.copy()
        o_plot = o.copy() * 1e6 / (r.copy() * BLASIUS_C)
        x_plot *= BLASIUS_C if useReynolds else 1.0
        axins.semilogy(x_plot, o_plot, linestyle=lin_fp, color=col_fp)

    if depthConstant:
        for r, o, x in zip(reynolds_bfs, omegas_bfs, x_bfs):
            x_plot = r.copy() if useReynolds else x.copy()
            o_plot = o.copy() * 1e6 / (r.copy() * BLASIUS_C)
            x_plot *= BLASIUS_C if useReynolds else 1.0
            axins.semilogy(x_plot, o_plot, linestyle=lin_bfs, color=col_bfs)

    for i, case in enumerate(cases.keys()):
        file_path = f"/home/victor/Desktop/PhD/data/neutralCurves/{case}/neutral_curve_temporal_gaster.dat"

        col = colors[i]
        # lin = linestyles[(i + 1) % len(linestyles)]
        lin = "-"

        reynolds, omegas, xs = read_neutral_curve(file_path)

        for r, o, x in zip(reynolds, omegas, xs):
            x_plot = r if useReynolds else x
            o_plot = o * 1e6 / (r * BLASIUS_C)
            x_plot *= BLASIUS_C if useReynolds else 1.0
            axins.semilogy(x_plot, o_plot, linestyle=lin, color=col, zorder=2)
    #
    # Zoom limits
    if useReynolds:
        axins.set_xlim(970, 1150)
    else:
        axins.set_xlim(re_to_x(970), re_to_x(1150))
    axins.set_ylim(40, 390)

    # Optional cosmetics
    axins.tick_params(labelsize=7)
    axins.set_ylim(40, 390)

    axins.tick_params(
        axis="y",
        which="major",
        left=True,
        right=False,
        labelleft=True,
        labelright=False,
    )

    axins.tick_params(
        axis="y",
        which="minor",
        left=True,
        right=False,
        labelleft=False,
        labelright=False,
    )

    # Draw connectors between inset and main plot
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.7", lw=0.5)

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.2),
        ncol=2,
    )
    # plt.show()
    plt.savefig(
        save_path,
        format="pdf",
        bbox_inches="tight",
    )


if __name__ == "__main__":
    main()
