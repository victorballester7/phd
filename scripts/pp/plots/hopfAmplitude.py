import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from pp.filterData import average_curves
from pp.hopfDataPoints import hopf2d_re_stable, hopf_re_unstable
from pp.figureFrame import FigureFrame, axis_on_top, jfm_ticks


plt.style.use("plots/style/jfm.mplstyle")

# gamma_1 and gamma_2 are read on the left axis (global norm), gamma_3 on the
# right one (norm over Omega'). The two groups are told apart by colour, and
# the right axis carries the colour of its curve: black and grey on the left,
# dark red on the right -- the same colours as the arrows of the stability
# panel.
COLORS = ["black", "0.45", "darkred"]
MARKERS = ["o", "s", "^"]
FIT_LS = ["--", "--", "--"]
COL3 = COLORS[2]

# gamma_3 amplitude on the right axis: right = K3_Y * left. The origin is the
# same on both axes, so the three paths share the bifurcation point.
K3_Y = 0.05


def mu_along(w, d, w_c, d_c, tangent):
    """
    Relative distance from the crossing point of a path,

        mu = +- sqrt((w/w_c - 1)^2 + (d/d_c - 1)^2),

    positive on the unstable side, so that the bifurcation is at mu = 0 on
    every path (eq. hopfdistance of the paper). Measuring w and d relative to
    their critical values puts the three paths on the same scale: along
    gamma_1 it is (w - w_c)/w_c, along gamma_3 (d - d_c)/d_c, whereas a
    distance in w units moved about thirty times faster along gamma_1 than
    one in d units along gamma_3.
    """
    side = np.sign((w - w_c) * tangent[0] + (d - d_c) * tangent[1])
    return side * np.hypot(w / w_c - 1.0, d / d_c - 1.0)


def makePlot():
    wd_c = np.array([[25.415, 3], [45.117, 1.736], [80, 1.405]])

    def bifLine2(w):
        return 0.0365 * w + 0.0898

    # gamma_1 and gamma_2: w, d, ||u'||_L2 and ||v'||_L2 over the whole domain
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
    ]

    # gamma_3 over the whole domain, same columns as above. Not plotted: past
    # the bifurcation the global norm is dominated by the wave the limit cycle
    # excites downstream of the gap, and it does not grow linearly (see text).
    gamma3_global = np.array(  # noqa: F841 -- kept for reference
        [
            [80, 1.41, 2.837676e-01, 1.787734e-01],
            [80, 1.42, 5.560265e-01, 3.447882e-01],
            [80, 1.43, 1.141631e00, 6.933020e-01],
            [80, 1.44, 6.076028e00, 3.787439e00],
            [80, 1.45, 8.950815e00, 5.609746e00],
            [80, 1.5, 1.194457e01, 7.804305e00],
            [80, 1.55, 1.398018e01, 8.999930e00],
        ]
    )

    # gamma_3 over Omega' = [w - 10, w] x [-d/2, 0.5], the end of the shear
    # layer over the gap up to the downstream edge. Columns: w, d,
    # sum over the history points in Omega' of A_u^2 + A_v^2, where A is half
    # the peak-to-peak amplitude over the last 15% of the run
    # (HistoryPoints.his of the d*_w80 DNS; five points: x = 70.6, 75.3 at
    # y = 0.5 and y = -d/2, and x = 80 at y = 0.5). The region stops at the
    # edge: the points downstream of it already carry the convected wave and
    # bend the curve upwards. Runs below d = 1.405 are still decaying towards
    # the equilibrium; beyond d = 1.45 the growth is no longer linear.
    gamma3_local = np.array(
        [
            [80, 1.35, 0],
            [80, 1.37, 0],
            [80, 1.38, 0],
            [80, 1.39, 0],
            [80, 1.40, 0],
            [80, 1.41, 6.011134e-03],
            [80, 1.42, 1.522502e-02],
            [80, 1.43, 2.435115e-02],
            [80, 1.44, 3.412315e-02],
            [80, 1.45, 4.607494e-02],
        ]
    )

    # ------------------------------------------------------------------
    # Figure with two panels
    # ------------------------------------------------------------------

    frame = FigureFrame(
        frame_w=4.6, aspect_ratio=0.85, pad_l=4.4, pad_b=0.85, pad_t=0.15
    )
    fig, ax = frame.fig, frame.ax

    # right panel: a small square in pad_r, vertically centred on the frame
    STAB_SIZE = 2.0  # cm
    # between the frame and the panel: the gamma_3 axis on the right of the
    # frame, then the panel's own d label and ticks
    STAB_GAP = 2.15  # cm
    HEIGHT_FROM_BOTTOM = 0.85
    ax_stab = frame.add_axes_cm(
        frame.pad_l + frame.frame_w + STAB_GAP,
        HEIGHT_FROM_BOTTOM,
        STAB_SIZE,
        STAB_SIZE,
    )

    # ------------------------------------------------------------------
    # LEFT PANEL : Hopf amplitude
    # ------------------------------------------------------------------

    bif_paths_tangent = [
        np.array([1.0, 0.0]),
        np.array([1.0, 0.0365]) / np.sqrt(1 + 0.0365**2),
        np.array([0.0, 1.0]),
    ]

    # (mu, A^2 in the units of its own axis, scale to the left axis)
    paths = []
    for i, block in enumerate(data):
        w_c, d_c = wd_c[i]
        mu = mu_along(block[:, 0], block[:, 1], w_c, d_c, bif_paths_tangent[i])
        paths.append((mu, block[:, 2] ** 2 + block[:, 3] ** 2, 1.0))
    w_c, d_c = wd_c[2]
    mu3 = mu_along(
        gamma3_local[:, 0], gamma3_local[:, 1], w_c, d_c, bif_paths_tangent[2]
    )
    paths.append((mu3, gamma3_local[:, 2], K3_Y))

    for i, (mu, amp2, scale) in enumerate(paths):
        c, m = COLORS[i], MARKERS[i]
        above = mu > 0.0
        fit_line = np.poly1d(np.polyfit(mu[above], amp2[above], 1))
        label = rf"$\boldsymbol\gamma_{{{i + 1}}}$"
        if i < 2:
            label += " (left axis)"
        if i == 2:
            label += " (right axis)"

        # Open markers on the stable side of the bifurcation: those points are
        # not small-amplitude limit cycles, they are equilibria with no limit
        # cycle at all, and a filled marker made them look like data lying on
        # a flat branch.
        ax.plot(mu[above], amp2[above] / scale, m, alpha=0.8, color=c, label=label)
        ax.plot(
            mu[~above],
            amp2[~above] / scale,
            m,
            alpha=0.8,
            markerfacecolor="none",
            markeredgecolor=c,
            markeredgewidth=0.7,
            linestyle="none",
        )

        # the fit is drawn until it leaves the frame
        x_fit = np.linspace(-0.005, 0.2, 200)
        y_fit = fit_line(x_fit) / scale
        keep = y_fit <= 1.6
        ax.plot(x_fit[keep], y_fit[keep], FIT_LS[i], color=c)

    ax.set_xlim([-0.08, 0.08])
    ax.set_ylim([-0.1, 1.6])
    ax.set_xticks([-0.05, 0, 0.05])
    ax.set_xlabel(r"$\mu$")

    ax.set_ylabel(
        r"$\||\boldsymbol{u}'|\|_{L_2(\Omega)}^{2}$",
        rotation=0,
        labelpad=10,
        va="center",  # horizontal labels on both sides: centre them on the axis
    )

    # inside the frame: the right of the plot holds the stability panel
    ax.legend(
        loc="center left", bbox_to_anchor=(1.33, 0.8), handlelength=1.0, frameon=False
        # loc="center left", bbox_to_anchor=(-0.85, 0.5), handlelength=1.0, frameon=False
    )
    ax.grid(False)
    frame.axis_on_top()
    frame.jfm_ticks()
    # the right edge belongs to the gamma_3 axis
    ax.tick_params(which="both", right=False)

    ax3 = ax.secondary_yaxis(
        "right", functions=(lambda y: y * K3_Y, lambda y: y / K3_Y)
    )
    ax3.set_yticks([0, 0.025, 0.050, 0.075])
    ax3.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax3.set_ylabel(
        r"$\||\boldsymbol{u}'|\|_{L_2(\Omega')}^{2}$",
        color=COL3,
        labelpad=4,
        rotation=0,
        va="center",
    )
    ax3.tick_params(which="both", colors=COL3)
    ax3.tick_params(which="minor", length=2)
    ax3.spines["right"].set_color(COL3)
    ax.spines["right"].set_color(COL3)  # drawn over the secondary spine

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

    alpha_fb = 0.2

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
    arrow_start = [(13, 3), (38, bifLine2(37)), (80, 1.0)]
    arrow_direction = [(1, 0), (1, 0.0365), (0, 1)]
    length = [31, 20, 0.95]
    labels = [
        r"$\boldsymbol\gamma_1$",
        r"$\boldsymbol\gamma_2$",
        r"$\boldsymbol\gamma_3$",
    ]

    for ars, ard, l, label, c in zip(
        arrow_start, arrow_direction, length, labels, COLORS
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

    ax_stab.set_xlim([0, 130])
    ax_stab.set_ylim([1, 4])

    ax_stab.set_xlabel(r"$w$", labelpad=2)
    ax_stab.set_ylabel(r"$d$", rotation=0, labelpad=7)

    # cleaner small panel
    ax_stab.set_xticks([0, 40, 80, 120])
    ax_stab.set_yticks([0, 1, 2, 3, 4])

    ax_stab.grid(False)
    axis_on_top(ax_stab)
    jfm_ticks(ax_stab)

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------

    script_path = os.path.dirname(os.path.abspath(__file__))

    save_path = os.path.join(
        script_path,
        "../../../images/hopfAmplitude.pdf",
    )

    fig.savefig(save_path, format="pdf")

    print("✓ Plot saved to:", save_path)


if __name__ == "__main__":
    makePlot()
