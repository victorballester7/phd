import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from pp.fft import fftFreqs
from pp.filterData import timeFilter
from pp.inputargs import parseArgs
from pp.colors import colors
from pp.fileManagement import extract_width_depth, readDataHistoryPointsMultiple
from pp.returnMap import getReturnPoints
from pp.analyzeData import compute_growthRate
from scipy.signal import hilbert
from scipy.stats import linregress
from scipy.fft import fft, fftfreq
from scipy.signal import butter, filtfilt
from scipy.signal import savgol_filter
from scipy.signal import find_peaks, savgol_filter
from scipy.interpolate import CubicSpline
from scipy.stats import linregress


def createFigure(returnMap: bool, is3d: bool):
    """
    Create a figure for plotting history points.

    Args:
        returnMap:(bool): Flag to indicate if the plot is for a return map.

    Returns:
        tuple: Figure and axes objects.
    """

    if returnMap:
        fig, axes = plt.subplots(1, 2, figsize=(8, 6))
        # axes = np.array([axes])  # Ensure axes is iterable
        # axes[0].set_aspect('equal', adjustable='box')
    else:
        numplots = 3 if not is3d else 4
        fig, axes = plt.subplots(1, numplots, figsize=(15, 5), sharex=True)

    for ax in axes:
        ax.grid(True)

    return fig, axes


def main():
    """
    Main function to plot history points from one or multiple folders as a function of time or in phase space (u vs v). The points are stored in the historyPoints.dat file in each folder.
    """
    args = parseArgs()

    folders = args.folders
    is3d = any("3d" in f.lower() for f in folders)

    data = readDataHistoryPointsMultiple(folders)

    fig, axes = createFigure(args.returnMap, is3d)

    onePointPerFolder = False
    if len(folders) > 1 and len(args.points) == len(folders):
        print(
            colors.WARNING
            + "Warning: Detected same number of folders and points. Would you like to plot each point from each folder? (Y/n)"
            + colors.ENDC
        )
        user_input = input().strip().lower()
        if user_input == "y" or user_input == "":
            onePointPerFolder = True
            print(colors.OKGREEN + "Plotting one point per folder." + colors.ENDC)
        else:
            print(
                colors.OKGREEN + "Plotting all points from all folders." + colors.ENDC
            )

    for j, f in enumerate(folders):
        isIncNS = "inc" in f.lower()
        isIncNSlabel = "(incNS)" if isIncNS else "(comNS)"
        is3d_f = "3d" in f.lower()
        points, time, fields = data[f]
        # fileds = np.swapaxes(fields, 0, 1)  # time, points, vars
        time, fields = timeFilter(time, fields, args.time_min, args.time_max)
        # fields = np.swapaxes(fields, 0, 1)  # points, time, vars

        depth, width = extract_width_depth(f)

        print(f"Processing folder: {f}")
        for i, p in enumerate(points):
            print(f"Point {i}: x = {p[0]}, y = {p[1]}, z = {p[2]}")

        print(colors.OKGREEN + "You are seeing the following points:" + colors.ENDC)
        if onePointPerFolder:
            p = args.points[j]
            print(
                colors.OKGREEN
                + f"x = {points[p, 0]}, y = {points[p, 1]}, z = {points[p, 2]} {isIncNSlabel}"
                + colors.ENDC
            )
        else:
            for p in args.points:
                print(
                    colors.OKGREEN
                    + f"x = {points[p, 0]}, y = {points[p, 1]}, z = {points[p, 2]} {isIncNSlabel}"
                    + colors.ENDC
                )
        print("\n")  # needed when calling with & disown

        for i, p in enumerate(args.points):
            if onePointPerFolder and i != j:
                continue  # Skip this point if we're plotting one point per folder and this isn't the corresponding point for this folder
            label = f"d{depth}_w{width}_x{points[p, 0]}_y{points[p, 1]}_z{points[p, 2]}"
            if not is3d_f:
                u = fields[p, :, 0] if isIncNS else fields[p, :, 1] / fields[p, :, 0]
                v = fields[p, :, 1] if isIncNS else fields[p, :, 2] / fields[p, :, 0]
                prho = fields[p, :, 2] if isIncNS else fields[p, :, 0]
                vars = (u, v, prho)
                labels_vars = ("u", "v", "p" if isIncNS else "ρ")
            else:
                u = fields[p, :, 0] if isIncNS else fields[p, :, 1] / fields[p, :, 0]
                v = fields[p, :, 1] if isIncNS else fields[p, :, 2] / fields[p, :, 0]
                w = fields[p, :, 2] if isIncNS else fields[p, :, 3] / fields[p, :, 0]
                prho = fields[p, :, 3] if isIncNS else fields[p, :, 0]
                vars = (u, v, w, prho)
                labels_vars = ("u", "v", "w", "p" if isIncNS else "ρ")
            # vars = (u, v, prho) if not is3d_f else (u, v, w, prho)
            if args.returnMap:
                var_section = u
                var_return = v if var_section is u else u
                section = 0.43
                min_section_return_variable = -0.3
                max_section_return_variable = -0.15
                # section_u = 0.44
                # min_section_v = -0.04
                # max_section_v = -0.01
                _, return_points = getReturnPoints(
                    time,
                    var_section,
                    var_return,
                    section,
                    min_section_return_variable,
                    max_section_return_variable,
                )  # return the values of time and the values of v corresponding to the return map with the section of u = const
                axes[0].plot(
                    u,
                    v,
                    label=f"u vs v (x = {points[p, 0]}, y = {points[p, 1]}, z = {points[p, 2]})",
                )
                axes[0].set_xlabel("u")
                axes[0].set_ylabel("v")

                # axes[1].plot(v_return[1:], v_return[:-1], 'o', label=f"Return map (x = {points[p, 0]}, y = {points[p, 1]})")
                n_points = len(return_points) - 1 if len(return_points) > 1 else 1
                scatter = True
                lab = f"Return map (x = {points[p, 0]}, y = {points[p, 1]}, z = {points[p, 2]})"
                if scatter:
                    # Scatter plot with color gradient
                    sm = plt.cm.ScalarMappable(
                        cmap="winter", norm=plt.Normalize(vmin=0, vmax=n_points)
                    )
                    col = plt.cm.winter(np.linspace(0, 1, n_points))  # rainbow colormap
                    plt.colorbar(sm, ax=axes[1], label="Point index (order)")
                    axes[1].scatter(
                        return_points[:-1],
                        return_points[1:],
                        s=10,
                        c=col,
                        cmap="winter",
                        label=lab,
                        alpha=0.7,
                    )
                else:
                    axes[1].plot(
                        return_points[1:],
                        return_points[:-1],
                        ".-",
                        label=lab,
                        alpha=0.7,
                    )

                axes[1].set_xlabel("v_n return point")
                axes[1].set_ylabel("v_n+1 return point")
                axes[1].set_xlim([np.min(var_return), np.max(var_return)])
                axes[1].set_ylim([np.min(var_return), np.max(var_return)])
            else:
                for f, lab, ax in zip(vars, labels_vars, axes):
                    if args.log:
                        f_plot = np.log(f - np.mean(f))
                    else:
                        f_plot = f

                    # def peak_envelope(t, u, smooth=True):
                    #     # Remove mean
                    #     u = u - np.mean(u)

                    #     # Find positive peaks
                    #     peaks, _ = find_peaks(u)

                    #     # Extract peak times and values
                    #     t_peaks = t[peaks]
                    #     u_peaks = u[peaks]

                    #     # Cubic spline interpolation
                    #     spline = CubicSpline(t_peaks, u_peaks)
                    #     A = spline(t)

                    #     # Optional light smoothing
                    #     if smooth:
                    #         window = int(0.05 * len(A))  # 5% of signal length
                    #         if window % 2 == 0:
                    #             window += 1
                    #         A = savgol_filter(A, window, 3)

                    #     return A

                    # def fit_stuart_landau(t, A):
                    #     dt = t[1] - t[0]

                    #     # Smooth A slightly to reduce derivative noise
                    #     window = int(0.03 * len(A))
                    #     if window % 2 == 0:
                    #         window += 1
                    #     A_smooth = savgol_filter(A, window, 3)

                    #     # Compute derivative
                    #     dA_dt = np.gradient(A_smooth, dt)

                    #     # Avoid division by very small A (early noise region)
                    #     mask = A_smooth > 0.05 * np.max(A_smooth)

                    #     x = A_smooth[mask]**2
                    #     y = (dA_dt[mask] / A_smooth[mask])

                    #     # Linear regression
                    #     slope, intercept, r_value, _, _ = linregress(x, y)

                    #     sigma = intercept
                    #     l = -slope

                    #     print(f"sigma = {sigma}")
                    #     print(f"l = {l}")
                    #     print(f"R^2 = {r_value**2}")

                    #     # Plot linear fit
                    #     # plt.scatter(x, y, s=10, alpha=0.5)
                    #     # plt.plot(x, intercept + slope*x, 'r', linewidth=2)
                    #     # plt.xlabel("A^2")
                    #     # plt.ylabel("(1/A) dA/dt")
                    #     # plt.title("Stuart–Landau Linear Fit")
                    #     # plt.show()

                    #     return sigma, l

                    # # Usage
                    # A = peak_envelope(time, f)
                    # sigma, l = fit_stuart_landau(time, A)
                    # print(A)
                    # print(A.shape)
                    # print(time.shape)
                    # print(f.shape)
                    # ax.plot(time, A + np.mean(f))

                    ax.plot(time, f_plot, label=label)

                    if args.growthrate:
                        growth_rate, intercept = compute_growthRate(time, f, ax)
                        print(
                            f"Estimated growth rate σ for {lab} at point {label}: {growth_rate:.6f}"
                        )
                        ax.plot(
                            time,
                            np.exp(growth_rate * time + intercept) + np.median(f),
                            linestyle="--",
                            label=f"Exp. growth (σ={growth_rate:.6f})",
                        )

                        args.fft = True

                    # if args.logyaxis:
                    #     ax.set_yscale("log")
                    # fit the stuart-landau equation to the data to estimate growth rate and frequency
                    # A(t) = R(t) * exp(i * theta(t))
                    # dA/dt = (sigma + i * omega) * A - l_r * |A|^2 * A
                    # dR/dt = sigma * R - l_r * R^3
                    # dtheta/dt = omega - l_i * R^2


                    # R(t)^-2 = (l_r / sigma) + (R0^-2 - l_r / sigma) * exp(-2 * sigma * t)
                    # theta(t) = theta0 + (omega - l_i / l_r * sigma) * t + l_i / l_r * ln(R(t)/R0)

                    # f_centered = f - np.mean(f)

                    # Find dominant frequency
                    # fft_vals = fft(f_centered)
                    # freqs = fftfreq(len(f_centered), time[1] - time[0])
                    # dominant_freq = np.abs(freqs[np.argmax(np.abs(fft_vals[1:len(freqs)//2]))+1])
                    #
                    # # Bandpass filter around dominant frequency
                    # fs = 1.0 / (time[1] - time[0])
                    # lowcut = dominant_freq * 0.5
                    # highcut = dominant_freq * 2.0
                    #
                    # nyq = 0.5 * fs
                    # b, a = butter(4, [lowcut/nyq, highcut/nyq], btype='band')
                    # filtered = filtfilt(b, a, f_centered)

                    # Hilbert transform
                    # f_hilbert = hilbert(f - np.mean(f))
                    # R = np.abs(f_hilbert)
                    # theta = np.unwrap(np.angle(f_hilbert))

                    # ax.plot(time, R + np.mean(f), label=f"Envelope {label} {lab}")
                    # ax.plot(time, np.gradient(theta, time), label=f"Phase {label} {lab}")

                    # # estimate omega
                    # slope, intercept, *_ = linregress(time, theta)
                    # omega_estimate = slope
                    # print(
                    #     f"Estimated angular frequency ω for {lab} at point {label}: {omega_estimate:.6f}"
                    # )

                    # omega_inst = np.gradient(theta, time)
                    # omega = np.mean(omega_inst)

                    #  if args.growthrate:
                    #      from scipy.signal import hilbert
                    #
                    #      t = time - time[0]               # shift so t starts at 0
                    #      f_c = f - np.mean(f)             # remove offset
                    #
                    #      # data-driven initial guesses
                    #      analytic = hilbert(f_c)
                    #      env = np.abs(analytic)
                    #      phase = np.unwrap(np.angle(analytic))
                    #      # sigma from slope of log(envelope)
                    #      mask = env > 0.01 * np.max(env)
                    #      sigma0, _, *_ = linregress(t[mask], np.log(env[mask]))
                    #      # omega from slope of unwrapped phase
                    #      omega0, _, *_ = linregress(t, phase)
                    #      a0_0 = env[0]
                    #      phi0 = phase[0]

                    #      print(f"Initial guesses for {lab} at point {label}: σ={sigma0:.6f}, ω={omega0:.6f}, φ={phi0:.6f}, A={a0_0:.6f}")
                    #
                    #      # plot envelope
                    #      ax.plot(time, np.log(env + np.mean(f)), linestyle="--", label=f"Envelope (Hilbert) {label} {lab}")
                    #
                    #
                    #      def function_to_fit(t, sigma, omega, phi, a0, b0):
                    #          return a0 * np.exp(sigma * t) * np.cos(omega * t + phi) + b0

                    #      # ax.plot(time, function_to_fit(t, sigma0, omega0, phi0, a0_0, np.mean(f)), linestyle="--", label=f"Initial guess fit: {a0_0:.3f}·exp({sigma0:.4f}·t)·cos({omega0:.4f}·t+{phi0:.2f})")
                    #      initial_guess = (sigma0, abs(omega0), phi0, a0_0, np.mean(f))
                    # try:
                    #     params_opt, _ = curve_fit(
                    #         function_to_fit, t, f,
                    #         p0=initial_guess,
                    #         bounds=([-0.1, 0, -np.pi, 0, -np.inf],
                    #                [0.1, 1, np.pi, np.inf, np.inf]),
                    #         maxfev=10000
                    #     )
                    #     sigma_opt, omega_opt, phi_opt, a0_opt, b0_opt = params_opt
                    #     print(f"σ={sigma_opt:.6f}, ω={omega_opt:.6f}, φ={phi_opt:.6f}, A={a0_opt:.6f}")
                    #     ax.plot(time, function_to_fit(t, *params_opt), linestyle=":",
                    #            label=f"Fit: {a0_opt:.3f}·exp({sigma_opt:.4f}·t)·cos({omega_opt:.4f}·t+{phi_opt:.2f})")
                    # except Exception as e:
                    #     print(f"Could not fit {lab} at {label}: {e}")

                    # # estimate growth rate
                    # slope_amp, intercept_amp, *_ = linregress(time, np.log(R))
                    # growth_rate_estimate = slope_amp
                    # print(
                    #     f"Estimated growth rate σ for {lab} at point {label}: {growth_rate_estimate:.6f}"
                    # )

                    # dRdt = np.gradient(R, time)
                    # y = dRdt / R
                    # x = R**2
                    # slope_gr, intercept_gr, *_ = linregress(x, y)
                    # sigma = intercept_gr
                    # l_r = -slope_gr
                    # print(
                    #     f"Estimated growth rate σ (from dR/dt) for {lab} at point {label}: {sigma:.6f}, non-linear coefficient l_r: {l_r:.6f}"
                    # )
                    # ax.plot(
                    #     time,
                    #     R,
                    #     label=f"Envelope {label} {lab}",
                    # )
                    # ax.plot(
                    #     time,
                    #     theta,
                    #     label=f"Phase {label} {lab}",
                    # )

                    ax.set_xlabel("t")
                    ax.set_ylabel(lab)

                    if args.fft:
                        time_new = time
                        f_new = f
                        if not isIncNS:
                            # timestep is not constant, cannot do proper fft
                            dt_array = np.diff(time)
                            dt_mean = np.mean(dt_array)
                            time_new = np.arange(time[0], time[-1], dt_mean)
                            f_new = np.interp(time_new, time, f)

                            ax.plot(
                                time_new,
                                f_new,
                                linestyle="--",
                                label=f"Interpolated {lab} for FFT",
                            )
                        # else:
                        #     # remove previous plot ax.plot
                        #     plt.close(fig)
                        fundamental_freqs = fftFreqs(
                            time_new, f_new, width, depth, points[p], lab
                        )

                        # if args.growthrate:
                        #     for omega0 in fundamental_freqs:
                        #         # multiply the signal by exp(-i * omega0 * t) to get the complex amplitude of the mode
                        #         f_complex = f * np.exp(-1j * omega0 * time)
                        #         fundamental_freqs = fftFreqs(time_new, f_complex, width, depth, points[p], lab)

                    elif args.fit and lab in ("u", "v", "w"):
                        signal = f

                        # fit to a * exp(growth_rate * t) * cos(omega * t + b)
                        def fit_func(t, a, growth_rate, omega, b):
                            return a * np.exp(growth_rate * t) * np.cos(omega * t + b)

                        # initial guess
                        params_init = [0.01, 0, 0.2, 0]  # a, growth_rate, omega, b
                        try:
                            params_opt, params_cov = curve_fit(
                                fit_func, time, signal, p0=params_init
                            )
                            a_opt, growth_rate_opt, omega_opt, b_opt = params_opt
                            print(
                                f"Fitted parameters for {lab} at point {label}: (A exp(σ t) cos(ω t + φ)"
                            )
                            print(
                                f"  A = {a_opt:.6f}, σ = {growth_rate_opt:.6f}, ω = {omega_opt:.6f}, φ = {b_opt:.6f}"
                            )
                            signal_fit = fit_func(time, *params_opt)

                            ax.plot(
                                time,
                                signal_fit,
                                linestyle=":",
                                label=f"Fit {label} {a_opt:.3f} exp({growth_rate_opt:.4f} t) cos({omega_opt:.4f} t + {b_opt:.2f})",
                            )
                        except Exception as _:
                            print(f"Could not fit data for {lab} at point {label}")

                    # growth_rate = 0.012835
                    # freq = 0.2190548
                    # a = 0.003
                    # phi = -1.3
                    # growth_rate = 0.0059853
                    # freq = 0.26979
                    # a = 0.002
                    # phi = -1.6
                    # growth_rate = 0.0047955
                    # freq = 0.0
                    # a = 0.8
                    # phi = 0
                    # ax.plot(
                    #     time,
                    #     growth_rate * time - 5,
                    #     linestyle=":",
                    #     label=f"Exp. growth (a={a}, σ={growth_rate}, f={freq})"
                    # )
                    # growth_rate = np.array(
                    #     [
                    #         0.000636,
                    #         0.0012473,
                    #         0.0020582,
                    #         0.0024620,
                    #         0.0026258,
                    #         0.0027111,
                    #         0.0027672,
                    #     ]
                    # )
                    # a = np.arange(len(growth_rate)) + 2
                    # labels = [
                    #     "T = 100, growth rate = 0.000636",
                    #     "T = 150, growth rate = 0.0012473",
                    #     "T = 300, growth rate = 0.0020582",
                    #     "T = 500, growth rate = 0.0024620",
                    #     "T = 700, growth rate = 0.0026258",
                    #     "T = 900, growth rate = 0.0027111",
                    #     "T = 1100, growth rate = 0.0027672",
                    # ]
                    # f = np.array(
                    #     [growth_rate[i] * time - a[i] for i in range(len(growth_rate))]
                    # ).T
                    # # plot all growth rates
                    # ax.plot(time, f, linestyle=":", label=labels)

    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
