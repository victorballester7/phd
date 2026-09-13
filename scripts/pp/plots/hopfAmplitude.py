import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from pp.fontSizeLaTex import compute_mpl_fontsize, compute_figure_size
from pp.colors import mix_with_black, colors
from pp.analyzeData import shared_depth_curve, build_bifurcation_pairs
from pp.smoothing import smooth_curve
from pp.filterData import average_curves
from pp.hopfDataPoints import hopf2d_re_stable, hopf_re_unstable

import matplotlib as mpl

plt.style.use("plots/style/jfm.mplstyle")


def makePlot():
    useVL2 = True

    wd_c = np.array([[25.415, 3], [45.117, 1.736], [80, 1.405]])

    uvL2_muc = np.array(
        [
            [407.9995, 0.4677355],
            [411.5315, 0.451751],
            [417.6885, 0.680634],
        ]
    )

    def bifLine2(w):
        return 0.0365 * w + 0.0898

    data = [
        np.array(
            [
                [24, 3, 0, 0],
                [25, 3, 0, 0],
                [25.415, 3, 0, 0],
                [25.55, 3, 1.911525e-01, 1.610369e-01],
                [25.6, 3, 2.125251e-01, 1.766622e-01],
                [25.8, 3, 2.954721e-01, 2.443431e-01],
                [26, 3, 3.648846e-01, 2.943190e-01],
                [27, 3, 5.669644e-01, 4.439065e-01],
                [30, 3, 1.001902e00, 7.250598e-01],
            ]
        ),
        np.array(
            [
                [43, bifLine2(43), 0, 0],
                [44, bifLine2(44), 0, 0],
                [45, 1.7316, 0, 0],
                [45.22, 1.74033, 1.059243e-01, 8.267704e-02],
                [45.5, 1.75055, 2.269154e-01, 1.729006e-01],
                [45.77, 1.7604, 2.930313e-01, 2.204655e-01],
                [46.2098, 1.7765, 3.550205e-01, 2.640201e-01],
                [46.3200, 1.7805, 3.846617e-01, 2.850623e-01],
                [46.4602, 1.7856, 4.131914e-01, 3.053044e-01],
                [46.8800, 1.8009, 4.538820e-01, 3.338911e-01],
                [47.1000, 1.8090, 4.802795e-01, 3.529278e-01],
            ]
        ),
        np.array(
            [
                [80, 1.35, 0, 0],
                [80, 1.37, 0, 0],
                [80, 1.38, 0, 0],
                [80, 1.39, 0, 0],
                [80, 1.4, 0, 0],
                [80, 1.41, 2.837676e-01, 1.787734e-01],
                [80, 1.42, 5.560265e-01, 3.447882e-01],
                [80, 1.43, 1.141631e00, 6.933020e-01],
                [80, 1.44, 6.076028e00, 3.787439e00],
                [80, 1.45, 8.950815e00, 5.609746e00],
                [80, 1.5, 1.194457e01, 7.804305e00],
                [80, 1.55, 1.398018e01, 8.999930e00],
            ]
        ),
    ]

    # ------------------------------------------------------------------
    # Figure with two panels
    # ------------------------------------------------------------------

    latex_width_cm = 10.0
    fig_width_in = latex_width_cm / 2.54

    fig = plt.figure(
        figsize=(fig_width_in, fig_width_in * 0.55),
        constrained_layout=False,
    )

    gs = GridSpec(
        3,
        2,
        width_ratios=[2, 1],
        height_ratios=[1, 2, 1],
        wspace=0.3,
    )

    ax = fig.add_subplot(gs[:, 0])  # left plot spans all rows
    ax_stab = fig.add_subplot(gs[1, 1])  # right plot only middle row

    # ------------------------------------------------------------------
    # LEFT PANEL : Hopf amplitude
    # ------------------------------------------------------------------

    bif_paths_tangent = [
        np.array([1.0, 0.0]),
        np.array([1.0, 0.0365]) / np.sqrt(1 + 0.0365**2),
        np.array([0.0, 1.0]),
    ]

    for i, (block, tangent) in enumerate(zip(data, bif_paths_tangent)):
        if i > 1:
            continue

        w, d = block[:, 0], block[:, 1]
        w_c, d_c = wd_c[i]

        mu = (w - w_c) * tangent[0] + (d - d_c) * tangent[1]

        idx = 2 if not useVL2 else 3

        idx_mu_positive = np.where(mu > 0)[0]

        mu_pos = mu[idx_mu_positive]
        block_pos = block[idx_mu_positive, idx]

        coeffs = np.polyfit(
            mu_pos,
            block_pos**2,
            1,
        )

        fit_line = np.poly1d(coeffs)

        x_fit = np.linspace(
            -0.2,
            max(3, np.max(mu_pos)),
            100,
        )

        col = plt.get_cmap("tab10")(i)

        ax.plot(
            mu,
            block[:, idx] ** 2,
            "o",
            alpha=0.7,
            color=col,
            label=rf"$\boldsymbol\gamma_{{{i + 1}}}$",
        )

        ax.plot(
            x_fit,
            fit_line(x_fit),
            "--",
            color=col,
        )

    ax.set_xlabel(r"$\mu-\mu_c$")

    ax.set_ylabel(
        r"$\|v'(\mu)\|_{L_2}^{2}$",
        # r"$\displaystyle\frac{\|v'(\mu)\|_{L_2}^{2}}{\|v'(\mu_c)\|_{L_2}^{2}}$",
        rotation=0,
        labelpad=20,
    )

    ax.legend()
    ax.grid(True)

    # ------------------------------------------------------------------
    # RIGHT PANEL : stability diagram
    # ------------------------------------------------------------------

    _, _, hopf2d = average_curves(
        hopf2d_re_stable[1000],
        hopf_re_unstable[1000],
        npts=500,
    )

    w_hopf = hopf2d[:, 0]
    d_hopf = hopf2d[:, 1]

    unstable_color = (0.85, 0.30, 0.05)

    w_plot = np.concatenate([w_hopf, [133]])
    d_plot = np.concatenate([d_hopf, [d_hopf[-1]]])

    w_extended = np.concatenate([[0], w_hopf, [140]])
    d_extended = np.concatenate([[5], d_hopf, [d_hopf[-1]]])

    alpha_fb = 0.3

    ax_stab.fill_between(
        w_extended,
        0,
        d_extended,
        color="tab:blue",
        alpha=alpha_fb,
        zorder=1,
    )

    ax_stab.fill_between(
        w_extended,
        d_extended,
        4.1,
        color=unstable_color,
        alpha=alpha_fb,
        zorder=1,
    )

    ax_stab.plot(
        w_plot,
        d_plot,
        color="tab:green",
        lw=1.5,
        zorder=3,
    )

    # draw an arrow to indicate the bifurcation direction
    # arrow_start = [(15, 3), (38, bifLine2(38)), (80, 0.9)]
    # arrow_direction = [(1, 0), (1,0.0365), (0, 1)]
    # length = [25, 20, 1.2]
    # labels = [r"$\boldsymbol\gamma_1$", r"$\boldsymbol\gamma_2$", r"$\boldsymbol\gamma_3$"]
    arrow_start = [(15, 3), (38, bifLine2(38)), (80, 0.9)]
    arrow_direction = [(1, 0), (1,0.0365), (0, 1)]
    length = [25, 20, 1.2]
    labels = [r"$\boldsymbol\gamma_1$", r"$\boldsymbol\gamma_2$", r"$\boldsymbol\gamma_3$"]
    colors = ["black", "black", "darkred"]

    for ars, ard, l, label, c in zip(
        arrow_start,
        arrow_direction,
        length,
        labels,
        colors
    ):
        end = (
            ars[0] + l * ard[0],
            ars[1] + l * ard[1],
        )

        ax_stab.annotate(
            "",
            xy=end,
            xytext=ars,
            arrowprops=dict(
                arrowstyle="->",
                color=c,
                lw=1.5,
            ),
        )
        
        

        ax_stab.text(
            ars[0],
            ars[1],
            label,
            color=c,
            ha="center",
            va="top",
        )

    ax_stab.set_xlim([0, 100])
    ax_stab.set_ylim([1, 4])

    ax_stab.set_xlabel(r"$w/\delta^*$")
    ax_stab.set_ylabel(r"$d/\delta^*$", rotation=0, labelpad=10)

    # cleaner small panel
    ax_stab.set_xticks([0, 25, 50, 75,100])
    ax_stab.set_yticks([0, 1, 2, 3, 4])

    # no legend
    # no points

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------

    script_path = os.path.dirname(os.path.abspath(__file__))

    save_path = os.path.join(
        script_path,
        "../../../images/hopfAmplitude.pdf",
    )

    plt.savefig(
        save_path,
        format="pdf",
        bbox_inches="tight",
    )

    print("✓ Plot saved to:", save_path)


if __name__ == "__main__":
    makePlot()
