import matplotlib.pyplot as plt
import numpy as np
from pp.fft import fftFreqs
from pp.filterData import timeFilter
from pp.inputargs import parseArgs
from pp.colors import colors
from pp.fileManagement import extract_depth_width, readErrorHistoryMultiple
from pp.returnMap import getReturnPoints


def createFigure(returnMap: bool):
    """
    Create a figure for plotting history points.

    Args:
        dynamicalSystem (bool): Flag to indicate if the plot is for a dynamical system.

    Returns:
        tuple: Figure and axes objects.
    """

    if returnMap:
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))
        # axes = np.array([axes])  # Ensure axes is iterable
    else:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=True)

    for ax in axes:
        ax.grid(True)

    return fig, axes


def main():
    """
    Main function to plot history points from one or multiple folders as a function of time or in phase space (u vs v). The points are stored in the historyPoints.dat file in each folder.
    """
    args = parseArgs()

    folders = args.folders

    data = readErrorHistoryMultiple(folders)

    fig, axes = createFigure(args.returnMap)

    for f in folders:
        isIncNS = "inc" in f.lower()
        isIncNSlabel = "(incNS)" if isIncNS else "(comNS)"
        is3d_f = "3d" in f.lower()
        

        time, errors = data[f]

        time = time - time[0]  # reset time to start at 0

        time, errors = timeFilter(time, errors, args.time_min, args.time_max)
        depth, width = extract_depth_width(f)

        uL2, vL2 = errors[:, 0], errors[:, 3]

        diff_vL2 = np.gradient(vL2, time)
        
        areaRectangle = 150*(100 + width + 1000)
        areaGap = width * depth
        area_domain = areaRectangle + areaGap
        print(f"Area of the domain for d{depth}_w{width}: {area_domain}")

        if args.returnMap:
            # customize to your needs
            vL2_returnMap = vL2
            uL2_returnMap = uL2
            var_section = vL2_returnMap
            var_return = vL2_returnMap if var_section is vL2_returnMap else vL2_returnMap
            section = 12.43
            min_section_return_variable = 417.09
            max_section_return_variable = 417.11
            ###################################

            _, return_points = getReturnPoints(
                time,
                var_section,
                var_return,
                section,
                min_section_return_variable,
                max_section_return_variable,
            )  # return the values of time and the values of vL2 corresponding to the return map with the section np.median(uL2)
            print(return_points)
            vL2_returnMap = vL2_returnMap - np.average(vL2_returnMap)
            uL2_returnMap = uL2_returnMap - np.average(uL2_returnMap)
            axes[0].plot(vL2_returnMap, uL2_returnMap, "-", label=f"d{depth}_w{width}")
            axes[0].set_xlabel("L2 error v/area_domain")
            axes[0].set_ylabel("L2 error u/area_domain")
            # axes[0].set_ylabel("diff L2 error v")

            n_points = len(return_points) - 1 if len(return_points) > 1 else 1
            scatter = True
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
                    label=f"Return map d{depth}_w{width}",
                    alpha=0.7,
                )
            else:
                axes[1].plot(
                    return_points[1:],
                    return_points[:-1],
                    ".-",
                    label=f"Return map d{depth}_w{width}",
                )
            axes[1].set_xlabel("nth return point")
            axes[1].set_ylabel("n+1-th return point")
            xmin = min(np.min(var_return), min_section_return_variable)
            xmax = max(np.max(var_return), max_section_return_variable)
            axes[1].set_xlim(xmin, xmax)
            axes[1].set_ylim(xmin, xmax)
        else:
            totalL2 = np.sqrt(uL2**2 + vL2**2)
            vars = [uL2, vL2, totalL2]
            labels_vars = ["u", "v", "sqrt(u^2 + v^2)"]
            axx = axes.flatten()
            for f, lab, ax in zip(vars, labels_vars, axx):
                ax.plot(time, f, label=f"d{depth}_w{width}")
                ax.set_xlabel("t")
                ax.set_ylabel(f"L2 error {lab}")

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
                    fftFreqs(time_new, f_new, width, depth, np.empty(3), lab)
           
            print(
                colors.OKBLUE
                + f"Mean u L2 error for d{depth}_w{width} from t={time[0]} to t={time[-1]}: {np.mean(uL2):.6e}"
                + colors.ENDC
            )
            print(
                colors.OKBLUE
                + f"Mean v L2 error for d{depth}_w{width} from t={time[0]} to t={time[-1]}: {np.mean(vL2):.6e}"
                + colors.ENDC
            )
            print(
                colors.OKBLUE
                + f"Mean total L2 error for d{depth}_w{width} from t={time[0]} to t={time[-1]}: {np.mean(totalL2):.6e}"
                + colors.ENDC
            )


    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
