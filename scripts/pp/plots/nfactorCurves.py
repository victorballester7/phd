import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from pp.DeltaN_computation import computeNx
from pp.fileManagement import extract_depth_width

RE_REFERENCE = 1000.0
BLASIUS_C = 1.7207876573

plt.style.use("plots/style/jfm.mplstyle")

linestyles = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1))]

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
    for i, block in enumerate(blocks):
        if not block.strip():
            raise ValueError(f"Block {i} is empty in file {file_path}")

        arr = np.loadtxt(block.splitlines())
        arrays.append(arr)

        reynolds.append(arr[:, 0])
        omegas.append(arr[:, 1])

    return reynolds, omegas


def main():
    depthConstant = False

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
        save_path = os.path.join(script_path, "../../../images/nfactorCurvesd1.5.pdf")
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
        save_path = os.path.join(script_path, "../../../images/nfactorCurvesw40.pdf")

    latex_width_cm = 6.5
    fig_width_in = latex_width_cm / 2.54
    _, ax = plt.subplots(
        figsize=(fig_width_in, fig_width_in * 0.75),
        constrained_layout=False,
    )
    re = 1000
    ylabel = r"$n$"
    # ylabel = r"$\omega_r$"

    def x_to_re(x):
        return re * np.sqrt(1.0 + BLASIUS_C**2 * np.asarray(x) / re)

    # flat plate case
    file_path_fp = "/home/victor/Desktop/PhD/src/flatPlateRe1000inc/directLinearSolver/blowingSuction/data/pointsavg_n600.dat"
    x_fp, Nx_fp = computeNx(file_path_fp, False)
    res_fp = x_to_re(x_fp)
    lin_fp = linestyles[1]

    # correct last data point due to numerical issues with 2D extrapolation
    Nx_fp[-1] = 2.0 * Nx_fp[-2] - Nx_fp[-3]

    ax.plot(res_fp, Nx_fp, linestyle=lin_fp, color=col_fp, label="Flat plate")

    if depthConstant:
        # do bfs as well
        file_path_bfs = f"/home/victor/Desktop/PhD/src/bfsRe{re}inc/directLinearSolver/blowingSuction/d1.5/data/pointsavg_n600.dat"
        x_bfs, Nx_bfs = computeNx(file_path_bfs, False)
        res_bfs = x_to_re(x_bfs)
        lin_bfs = linestyles[2]
        # lin_bfs = "-"
        ax.plot(res_bfs, Nx_bfs, linestyle=lin_bfs, color=col_bfs, label="BFS")

    for i, case in enumerate(cases.keys()):
        file_path = f"/home/victor/Desktop/PhD/src/incGapRe{re}/directLinearSolver/blowingSuction/{case}/data/pointsavg_n600.dat"
        x, Nx = computeNx(file_path, False)
        res = x_to_re(x)
        col = colors[i]
        lin = linestyles[(i + 1) % len(linestyles)]
        lin = "-"
        ax.plot(res, Nx, linestyle=lin, color=col, label=rf"${cases[case]}\delta^*$")
        _, w = extract_depth_width(case)
        # filter x between w+80 and w+80+75
        x_start =  w + 80 
        x_end = x_start + 50
        x_avg_idx = (x >= x_start) & (x <= x_end)
        Nx_avg = Nx[x_avg_idx]
        re_avg = res[x_avg_idx]
        col_darken = col * 0.85
        ax.plot(
            re_avg,
            Nx_avg,
            linestyle="",
            marker="|",
            markersize=3,
            color=col_darken
        )

    ax.set_xlabel(r"$\mbox{\textit{Re}}$")
    ax.set_ylabel(ylabel, rotation=0, labelpad=10)

    def re_to_x(re):
        return (
            ((np.asarray(re) / RE_REFERENCE) ** 2 - 1.0) * RE_REFERENCE / BLASIUS_C**2
        )

    ax2 = ax.secondary_xaxis(
        "top",
        functions=(re_to_x, x_to_re),
    )
    ax2.set_xlabel(r"$x/\delta^*$")

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
