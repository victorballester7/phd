import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator, ScalarFormatter
from pp.DeltaN_computation import computeAmplitude, computeNx
from pp.fileManagement import extract_depth_width
from pp.figureFrame import FigureFrame
from pp.fileManagement import readFieldsBySection

RE_REFERENCE = 1000.0
BLASIUS_C = 1.7207876573

plt.style.use("plots/style/jfm.mplstyle")

linestyles = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1))]

# Increasing w gets darker. The reversed ramp used before made w = 80, the
# case the discussion is about, the palest curve on the plot.
cmap_tmp = plt.get_cmap("Purples")

DNS_2D = "/home/victor/Desktop/PhD/src/incGapRe1000/baseflow/dns"
DNS_3D = "/home/victor/Desktop/PhD/src/quasi3DGapRe1000inc/baseflow/dns"

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

# Amplitude at which the TS wave breaks down into a 3D flow, set by hand from
# the w = 66 case: the quasi-3D run leaves the 2D run (by 5%) at x - w ~ 385,
# where ||eta|| ~ 0.27. The w = 70 case leaves its 2D run at ~0.24, close
# enough for a single threshold.
ETA_BREAKDOWN = 0.27


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


def spanwise_norm(x, y, data):
    """L2 norm over y of the spanwise velocity w, the 3D-ness indicator.

    w is identically zero in a 2D flow, so any non-zero value marks the
    breakdown into a three-dimensional state. Not smoothed: its fluctuation
    along x is what the plot shows. Note that the w column of the
    points file is sampled on a single spanwise plane, unlike ``uvw``.
    """
    ymax = 150
    yindx = np.argmin(np.abs(y - ymax))
    return np.sqrt(np.trapezoid(data[:, :yindx, 4] ** 2, y[:yindx], axis=1))


def main():
    script_path = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(script_path, "../../../images/L2POvsX.pdf")
    # one curve per width, all on the same colour ramp: the 2D runs up to
    # w = 65, the quasi-3D runs for w = 70 and 80
    cases = {
        "d1.5_w60": (DNS_2D, "w=60"),
        "d1.5_w63": (DNS_2D, "w=63"),
        "d1.5_w64": (DNS_2D, "w=64"),
        "d1.5_w65": (DNS_2D, "w=65"),
        # "d1.5_w66": (DNS_2D, "w=66"),
        # "d1.5_w67": (DNS_2D, "w=67"),
        # "d1.5_w70": (DNS_2D, "w=70"),
        "d1.5_w66_Lz6.0_k32": (DNS_3D, "w=66"),
        "d1.5_w70_Lz6.0_k32": (DNS_3D, "w=70"),
        # "d1.5_w80_Lz6.0_k32": (DNS_3D, "w=80"),
    }
    # 2D runs drawn dashed under the quasi-3D run of the same width, in its
    # colour: the 2D oscillation the flow follows if it is kept
    # spanwise-invariant
    cases_2d_ref = {
        "d1.5_w66_Lz6.0_k32": "d1.5_w66",
        "d1.5_w70_Lz6.0_k32": "d1.5_w70",
    }

    n = len(cases)
    colors = cmap_tmp(np.linspace(0.35, 1.0, n))
    # the spanwise velocity is drawn behind the eta curves in a contrasting
    # colour: where the flow stays 2D it has zero width and only the purple
    # curve is visible. Teal, so that red is free for the breakdown amplitude
    col_w = "#1b9e77"
    col_breakdown = "tab:red"

    frame = FigureFrame(frame_w=5.6, aspect_ratio=0.75, pad_l=2.75, pad_b=0.8, pad_t=0.1)
    fig, ax = frame.fig, frame.ax
    ylabel = r"$\|\eta\|_{L^2(\mathcal{Y}(x))}$"
    # ylabel = r"$\omega_r$"

    for i, (case, (root, label)) in enumerate(cases.items()):
        file_path = f"{root}/{case}/data/pointsPO_n600.dat"
        if root == DNS_2D:
            # the 2D runs are a 3D flow with w = 0 everywhere
            x, L2 = computeAmplitude(file_path, False, field="|u|")
            L2_w = np.zeros_like(L2)
        else:
            x, y, data = readFieldsBySection(file_path)
            L2_w = spanwise_norm(x, y, data)
            L2 = np.sqrt(np.abs(data[:, :, 6]))

            ymax = 150
            yindx = np.argmin(np.abs(y - ymax))

            L2 = np.sqrt(np.trapezoid(L2[:, :yindx] ** 2, y[:yindx], axis=1))

        # window average L2 norm to smooth out the curve
        L2_smooth = gaussian_filter1d(L2, sigma=10)
        # L2_smooth = np.convolve(L2, np.ones(window_size) / window_size, mode="same")

        col = colors[i]
        # lin = linestyles[i]
        lin = "-"
        _, width = extract_depth_width(case)
        x = x - width
        L2, L2_smooth = suppress_unconverged_tail(case, x, L2, L2_smooth)

        # instantaneous ||w|| as a band around the averaged eta curve, in the
        # same units as eta so no second axis is needed; clipped below so the
        # lower edge stays on the log axis
        lower = np.maximum(L2_smooth - L2_w, 0.5 * L2_smooth)
        ax.fill_between(
            x, lower, L2_smooth + L2_w, color=col_w, alpha=0.55, linewidth=0, zorder=1
        )
        ax.plot(x, L2_smooth, linestyle=lin, color=col, label=rf"${label}$", zorder=2)

        if case in cases_2d_ref:
            file_path_2d = f"{DNS_2D}/{cases_2d_ref[case]}/data/pointsPO_n600.dat"
            x_2d, L2_2d = computeAmplitude(file_path_2d, False, field="|u|")
            L2_2d_smooth = gaussian_filter1d(L2_2d, sigma=10)
            ax.plot(x_2d - width, L2_2d_smooth, linestyle="--", color=col, zorder=1.5)

 #    # # for i, case in enumerate(cases3D.keys()):
 #    file_path = "/home/victor/Desktop/PhD/src/quasi3DGapRe1000inc/baseflow/dns/d1.5_w70_Lz6.0_k32/data/pointsPO_n600.dat"
 #    x, L2 = computeAmplitude(file_path, False, field="|u|")
 #    x, y, data = readFieldsBySection(file_path)
 #    L2 = np.sqrt(np.abs(data[:, :, 6]))
 #
 #    # integrate from 0 to ymax
 #    ymax = 150
 #    yindx = np.argmin(np.abs(y - ymax))
 #
 #    L2 = np.sqrt(np.trapezoid(L2[:, :yindx] ** 2, y[:yindx]))
 #
 #
 #    # window average L2 norm to smooth out the curve
 #    L2_smooth = gaussian_filter1d(L2, sigma=10)
 #    # L2_smooth = np.convolve(L2, np.ones(window_size) / window_size, mode="same")
 #
 #    col = "tab:red"
 #    # lin = linestyles[i]
 #    lin = "-"
 #    _, width = extract_depth_width("d1.5_w70")
 #    x = x - width
 #    L2, L2_smooth = suppress_unconverged_tail("d1.5_w70", x, L2, L2_smooth)
 #    ax.plot(x, L2, linestyle="-", color=col, linewidth=0.6, alpha=0.4)
 #    ax.plot(x, L2_smooth, linestyle=lin, color=col)
 #    ax.plot(x, L2_smooth, linestyle=lin, color=col, label=r"$w=70$") 
 #
 # # # for i, case in enumerate(cases3D.keys()):
 #    file_path = "/home/victor/Desktop/PhD/src/quasi3DGapRe1000inc/baseflow/dns/d1.5_w70_Lz6.0_k32/data/pointsPO_n600.dat"
 #    x, L2 = computeAmplitude(file_path, False, field="|u|")
 #    x, y, data = readFieldsBySection(file_path)
 #    L2 = np.sqrt(np.abs(data[:, :, 6]))
 #
 #    # integrate from 0 to ymax
 #    ymax = 150
 #    yindx = np.argmin(np.abs(y - ymax))
 #
 #    L2 = np.trapezoid(L2[:, :yindx] ** 2, y[:yindx])
 #
 #
 #    # window average L2 norm to smooth out the curve
 #    L2_smooth = gaussian_filter1d(L2, sigma=10)
 #
 #    L2_smooth = np.sqrt(L2_smooth)  # take the square root after smoothing the squared values
 #    # L2_smooth = np.convolve(L2, np.ones(window_size) / window_size, mode="same")
 #
 #    col = "tab:green"
 #    # lin = linestyles[i]
 #    lin = "-"
 #    _, width = extract_depth_width("d1.5_w70")
 #    x = x - width
 #    L2, L2_smooth = suppress_unconverged_tail("d1.5_w70", x, L2, L2_smooth)
 #    ax.plot(x, L2, linestyle="-", color=col, linewidth=0.6, alpha=0.4)
 #    ax.plot(x, L2_smooth, linestyle=lin, color=col)
 #    ax.plot(x, L2_smooth, linestyle=lin, color=col, label=r"$w=70$") 


    # file_path = "/home/victor/Desktop/PhD/src/quasi3DGapRe1000inc/baseflow/dns/d1.5_w65_Lz6.0_k32/data/pointsPO_n600.old.dat"
    # x, L2 = computeAmplitude(file_path, False, field="|u|")
    # # window average L2 norm to smooth out the curve
    # L2_smooth = gaussian_filter1d(L2, sigma=10)
    # # L2_smooth = np.convolve(L2, np.ones(window_size) / window_size, mode="same")
    #
    # col = "tab:green"
    # # lin = linestyles[i]
    # lin = "-"
    # _, width = extract_depth_width("d1.5_w65")
    # x = x - width
    # L2, L2_smooth = suppress_unconverged_tail("d1.5_w65", x, L2, L2_smooth)
    # ax.plot(x, L2, linestyle="-", color=col, linewidth=0.6, alpha=0.4)
    # ax.plot(x, L2_smooth, linestyle=lin, color=col)
    # ax.plot(x, L2_smooth, linestyle=lin, color=col, label=r"old$w=65$") 
    #
    # for i, case in enumerate(cases3D.keys()):
    #     file_path = f"/home/victor/Desktop/PhD/src/quasi3DGapRe1000inc/baseflow/dns/{case}/data/pointsPO_n600.oold.dat"
    #     x, L2 = computeAmplitude(file_path, False, field="|u|")
    #     # window average L2 norm to smooth out the curve
    #     L2_smooth = gaussian_filter1d(L2, sigma=10)
    #     # L2_smooth = np.convolve(L2, np.ones(window_size) / window_size, mode="same")
    #
    #     col = plt.get_cmap("Oranges")(np.linspace(0.45, 1.0, n))[i]
    #     # lin = linestyles[i]
    #     lin = "-"
    #     _, width = extract_depth_width(case)
    #     x = x - width
    #     L2, L2_smooth = suppress_unconverged_tail(case, x, L2, L2_smooth)
    #     for a in (ax, axins):
    #         a.plot(x, L2, linestyle="-", color=col, linewidth=0.6, alpha=0.4)
    #     axins.plot(x, L2_smooth, linestyle=lin, color=col)
    #     ax.plot(x, L2_smooth, linestyle=lin, color=col, label=rf"oold${cases3D[case]}$")
    #
    #
    # for i, case in enumerate(cases3D.keys()):
    #     file_path = f"/home/victor/Desktop/PhD/src/quasi3DGapRe1000inc/baseflow/dns/{case}/data/pointsPO_n600.ooold.dat"
    #     x, L2 = computeAmplitude(file_path, False, field="|u|")
    #     # window average L2 norm to smooth out the curve
    #     L2_smooth = gaussian_filter1d(L2, sigma=10)
    #     # L2_smooth = np.convolve(L2, np.ones(window_size) / window_size, mode="same")
    #
    #     col = plt.get_cmap("Blues")(np.linspace(0.45, 1.0, n))[i]
    #     # lin = linestyles[i]
    #     lin = "-"
    #     _, width = extract_depth_width(case)
    #     x = x - width
    #     L2, L2_smooth = suppress_unconverged_tail(case, x, L2, L2_smooth)
    #     for a in (ax, axins):
    #         a.plot(x, L2, linestyle="-", color=col, linewidth=0.6, alpha=0.4)
    #     axins.plot(x, L2_smooth, linestyle=lin, color=col)
    #     ax.plot(x, L2_smooth, linestyle=lin, color=col, label=rf"ooold${cases3D[case]}$")


    ax.set_xlabel(r"$x - w$")
    ax.set_xlim(0, 900)
    # logarithmic: the norm spans three decades, and on a linear axis the two
    # damped cases -- the ones the text describes as decaying to negligible
    # values -- were indistinguishable from a flat line at zero
    ax.set_yscale("log")
    ax.set_ylim(3e-4, 1.2)
    ax.set_ylabel(ylabel, rotation=0, labelpad=4, ha="right", va="center")

    ax.grid(False)
    frame.axis_on_top()
    frame.jfm_ticks()
    # no minor ticks on the log axis: one per unit inside each decade crowds
    # the three decades shown
    ax.yaxis.set_minor_locator(NullLocator())

    # legend entry for the dashed 2D curves, in a neutral colour since each one
    # takes the colour of its width
    ax.plot([], [], linestyle="--", color="0.4", label="2D")

    # breakdown amplitude across the whole frame, behind the eta curves
    ax.axhline(
        ETA_BREAKDOWN, color=col_breakdown, linestyle=":", linewidth=1.0,
        zorder=1.8, label="breakdown",
    )

    # legend entry for the band, after the eta curves: its half-width is the
    # wall-normal norm of zeta = <u_z'^2(t*)>_Z^(1/2), the spanwise rms of u_z'
    # at a single instant t*, defined in the caption as the counterpart of eta
    ax.fill_between(
        [], [], color=col_w, alpha=0.55, linewidth=0,
        label=r"$\|\zeta\|_{L^2(\mathcal{Y}(x))}$",
    )

    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
    )
    # plt.show()
    fig.savefig(save_path, format="pdf")


if __name__ == "__main__":
    main()
