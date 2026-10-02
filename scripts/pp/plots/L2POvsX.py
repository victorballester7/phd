import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator, ScalarFormatter
from matplotlib.transforms import TransformedBbox
from mpl_toolkits.axes_grid1.inset_locator import BboxConnector, BboxPatch
from pp.DeltaN_computation import computeAmplitude, computeNx
from pp.fileManagement import extract_depth_width
from pp.figureFrame import FigureFrame, axis_on_top

RE_REFERENCE = 1000.0
BLASIUS_C = 1.7207876573

plt.style.use("plots/style/jfm.mplstyle")

linestyles = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1))]

# Increasing w gets darker. The reversed ramp used before made w = 80, the
# case the discussion is about, the palest curve on the plot.
n = 6
cmap_tmp = plt.get_cmap("Purples")
colors = cmap_tmp(np.linspace(0.35, 1.0, n))
cmap_3d_tmp = plt.get_cmap("Reds")
colors_3d = cmap_3d_tmp(np.linspace(0.45, 1.0, n))

# ---------------------------------------------------------------------------
# Unconverged tail of the w = 60 run.
#
# That simulation has not converged past x - w ~ 650: the norm bottoms out at
# 6.5e-4 and then climbs by a factor of six to the end of the domain, which is
# a numerical artefact and not boundary-layer growth -- the frequency of this
# case lies above the upper branch everywhere downstream, so the wave can only
# decay. The tail is therefore replaced by a continuation of the exponential
# decay the run does resolve, fitted over the window below.
#
# The replacement is a model, not data. It has to be stated wherever the figure
# is used.
# ---------------------------------------------------------------------------
TAIL_FIX_CASE = "d1.5_w60"
TAIL_FIT_WINDOW = (450.0, 650.0)  # where the decay rate is measured
# in x - w: the run is handed over to the fitted decay across this window.
# A hard switch at one station left a step (the fit does not pass exactly
# through the data there), and scaling the instantaneous trace by the ratio
# to the unconverged mean carried that run's drift into the tail as wiggles.
# The window is centred on x - w ~ 610, where the run crosses the fit: the two
# agree to within 2% across it, while past 640 the run already drifts away
# (5% low by 700) and blending over that stretch left an S in the curve.
TAIL_BLEND = (580.0, 640.0)


def _smoothstep(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def suppress_unconverged_tail(case, x, L2, L2_smooth):
    """Hand the unconverged tail over to the fitted exponential decay.

    Both curves are blended in log space, with a smoothstep weight across
    ``TAIL_BLEND``: the smoothed curve goes from the run to the fitted decay,
    and the instantaneous trace goes from the run to that corrected smoothed
    curve. Downstream of the window both are the fitted decay.
    """
    if case != TAIL_FIX_CASE:
        return L2, L2_smooth

    fit = (x >= TAIL_FIT_WINDOW[0]) & (x <= TAIL_FIT_WINDOW[1])
    slope, intercept = np.polyfit(x[fit], np.log(L2_smooth[fit]), 1)
    log_decay = intercept + slope * x

    s = _smoothstep((x - TAIL_BLEND[0]) / (TAIL_BLEND[1] - TAIL_BLEND[0]))
    log_smooth = (1.0 - s) * np.log(L2_smooth) + s * log_decay
    log_inst = (1.0 - s) * np.log(L2) + s * log_smooth
    return np.exp(log_inst), np.exp(log_smooth)


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(script_path, "../../../images/L2POvsX.pdf")
    cases = {
        "d1.5_w60": "w=60",
        "d1.5_w63": "w=63",
        "d1.5_w64": "w=64",
        "d1.5_w65": "w=65",
        "d1.5_w70": "w=70",
        "d1.5_w80": "w=80",
    }

    cases3D = {
        "d1.5_w60_Lz5_k32": "w=60",
        "d1.5_w63_Lz6.0_k32": "w=63",
        "d1.5_w65_Lz6.0_k32": "w=65",
        "d1.5_w70_Lz6.0_k32": "w=70",
        "d1.5_w80_Lz6.0_k32": "w=80",
    }

    frame = FigureFrame(frame_w=5.6, aspect_ratio=0.75, pad_l=3.2, pad_b=0.8, pad_t=0.1)
    fig, ax = frame.fig, frame.ax
    re = 1000
    ylabel = r"$\||\boldsymbol{u}'|\|_{L^2(\mathcal{Y}(x))}$"
    # ylabel = r"$\omega_r$"

    # Zoom on the first hundred units behind the gap edge, where every case
    # still carries the shape of the oscillation shed by the downstream corner
    # and the curves have not yet separated.
    # The panel sits outside the frame, in the bottom of the left padding: its
    # bottom is level with the frame's and its top stays below the y label,
    # which is centred on the frame
    INS_W, INS_H = 1.8, 1.8  # cm
    INS_LEFT = 0.65  # cm from the left figure edge (room for its tick labels)
    axins = frame.add_axes_cm(INS_LEFT, frame.pad_b, INS_W, INS_H)

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
        L2, L2_smooth = suppress_unconverged_tail(case, x, L2, L2_smooth)
        for a in (ax, axins):
            a.plot(x, L2, linestyle="-", color=col, linewidth=0.6, alpha=0.4)
        axins.plot(x, L2_smooth, linestyle=lin, color=col)
        ax.plot(x, L2_smooth, linestyle=lin, color=col, label=rf"${cases[case]}$")

    for i, case in enumerate(cases3D.keys()):
        file_path = f"/home/victor/Desktop/PhD/src/quasi3DGapRe1000inc/baseflow/dns/{case}/data/pointsPO_n600.dat"
        x, L2 = computeAmplitude(file_path, False, field="|u|")
        # window average L2 norm to smooth out the curve
        L2_smooth = gaussian_filter1d(L2, sigma=10)
        # L2_smooth = np.convolve(L2, np.ones(window_size) / window_size, mode="same")

        col = colors_3d[i]
        # lin = linestyles[i]
        lin = "-"
        _, width = extract_depth_width(case)
        x = x - width
        L2, L2_smooth = suppress_unconverged_tail(case, x, L2, L2_smooth)
        for a in (ax, axins):
            a.plot(x, L2, linestyle="-", color=col, linewidth=0.6, alpha=0.4)
        axins.plot(x, L2_smooth, linestyle=lin, color=col)
        ax.plot(x, L2_smooth, linestyle=lin, color=col, label=rf"${cases3D[case]}$")
    
    for i, case in enumerate(cases3D.keys()):
        file_path = f"/home/victor/Desktop/PhD/src/quasi3DGapRe1000inc/baseflow/dns/{case}/data/pointsPO_n600.old.dat"
        x, L2 = computeAmplitude(file_path, False, field="|u|")
        # window average L2 norm to smooth out the curve
        L2_smooth = gaussian_filter1d(L2, sigma=10)
        # L2_smooth = np.convolve(L2, np.ones(window_size) / window_size, mode="same")

        col = plt.get_cmap("Greens")(np.linspace(0.45, 1.0, n))[i]
        # lin = linestyles[i]
        lin = "-"
        _, width = extract_depth_width(case)
        x = x - width
        L2, L2_smooth = suppress_unconverged_tail(case, x, L2, L2_smooth)
        for a in (ax, axins):
            a.plot(x, L2, linestyle="-", color=col, linewidth=0.6, alpha=0.4)
        axins.plot(x, L2_smooth, linestyle=lin, color=col)
        ax.plot(x, L2_smooth, linestyle=lin, color=col, label=rf"${cases3D[case]}$")


    for i, case in enumerate(cases3D.keys()):
        file_path = f"/home/victor/Desktop/PhD/src/quasi3DGapRe1000inc/baseflow/dns/{case}/data/pointsPO_n600.old.old.dat"
        x, L2 = computeAmplitude(file_path, False, field="|u|")
        # window average L2 norm to smooth out the curve
        L2_smooth = gaussian_filter1d(L2, sigma=10)
        # L2_smooth = np.convolve(L2, np.ones(window_size) / window_size, mode="same")

        col = plt.get_cmap("Blues")(np.linspace(0.45, 1.0, n))[i]
        # lin = linestyles[i]
        lin = "-"
        _, width = extract_depth_width(case)
        x = x - width
        L2, L2_smooth = suppress_unconverged_tail(case, x, L2, L2_smooth)
        for a in (ax, axins):
            a.plot(x, L2, linestyle="-", color=col, linewidth=0.6, alpha=0.4)
        axins.plot(x, L2_smooth, linestyle=lin, color=col)
        ax.plot(x, L2_smooth, linestyle=lin, color=col, label=rf"${cases3D[case]}$")


    ax.set_xlabel(r"$x - w$")
    ax.set_xlim(0, 900)
    # logarithmic: the norm spans three decades, and on a linear axis the two
    # damped cases -- the ones the text describes as decaying to negligible
    # values -- were indistinguishable from a flat line at zero
    ax.set_yscale("log")
    ax.set_ylim(3e-4, 1.2)
    ax.set_ylabel(ylabel, rotation=0, labelpad=30)

    axins.set_xlim(0, 75)
    axins.set_ylim(0.01, 0.18)

    # zoomed region on the main axes, joined to the panel on its left:
    # mark_inset can only join matching corners, so the connectors are built
    # by hand (panel upper/lower right -> region upper/lower left)
    zoom_kw = dict(fc="none", ec="0.7", lw=0.5)
    zoom = TransformedBbox(axins.viewLim, ax.transData)
    ax.add_patch(BboxPatch(zoom, **zoom_kw))
    # the connectors belong to ax, below its axis artists (zorder 10 from
    # axis_on_top), so they pass under the y tick labels and the y label
    for loc_ins, loc_zoom in ((2, 2), (4, 4)):
        connector = BboxConnector(
            axins.bbox, zoom, loc1=loc_ins, loc2=loc_zoom, zorder=5, **zoom_kw
        )
        connector.set_clip_on(False)
        ax.add_patch(connector)
    axins.tick_params(labelsize=7)
    axins.grid(False)
    axins.patch.set_alpha(1.0)

    ax.grid(False)
    frame.axis_on_top()
    axis_on_top(axins)
    frame.jfm_ticks()
    # the inset keeps its major ticks, mirrored, but no minor ones
    axins.tick_params(top=True, right=True)
    # no minor ticks on the log axis: one per unit inside each decade crowds
    # the three decades shown
    ax.yaxis.set_minor_locator(NullLocator())

    # white behind the y labels, so the connectors read as passing beneath
    # them rather than through them
    label_bg_ax = dict(facecolor="white", edgecolor="none", pad=0.5)
    label_bg_tick = dict(facecolor="white", edgecolor="none", pad=3)
    ax.yaxis.label.set_bbox(label_bg_ax)
    for label in ax.get_yticklabels():
        label.set_bbox(label_bg_tick)

    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
    )
    plt.show()
    fig.savefig(save_path, format="pdf")


if __name__ == "__main__":
    main()
