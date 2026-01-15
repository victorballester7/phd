import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from pp.fft import fftFreqs
from pp.filterData import timeFilter
from pp.inputargs import parseArgs
from pp.colors import colors
from pp.fileManagement import extract_width_depth, readDataHistoryPointsMultiple
from pp.returnMap import getReturnPoints
from scipy.signal import hilbert


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

    for f in folders:
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
        for p in args.points:
            print(
                colors.OKGREEN
                + f"x = {points[p, 0]}, y = {points[p, 1]}, z = {points[p, 2]} {isIncNSlabel}"
                + colors.ENDC
            )
        print("\n")  # needed when calling with & disown

        for i, p in enumerate(args.points):
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
                    ax.plot(time, f_plot, label=label)
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
                            
                            ax.plot(time_new, f_new, linestyle="--", label=f"Interpolated {lab} for FFT")
                        else:
                            # remove previous plot ax.plot
                            plt.close(fig) 
                        fftFreqs(time_new, f_new, width, depth, points[p], lab)
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
                            print(
                                f"Could not fit data for {lab} at point {label}"
                            )

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
