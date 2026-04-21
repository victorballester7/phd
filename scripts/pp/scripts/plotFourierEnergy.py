import matplotlib.pyplot as plt
import numpy as np
from pp.fft import fftFreqs
from pp.filterData import timeFilter
from pp.inputargs import parseArgs
from pp.colors import colors
from pp.fileManagement import (
    extract_width_depth,
    readFourierEnergyMultiple,
    extractValueXML,
)
from pp.returnMap import getReturnPoints
import os

# mix tab20b and tab20c to get more colors
tableau_colors = plt.get_cmap("tab20b").colors + plt.get_cmap("tab20c").colors
line_styles = ["-", "--", "-.", ":"]


def main():
    """
    Main function to plot history points from one or multiple folders as a function of time or in phase space (u vs v). The points are stored in the historyPoints.dat file in each folder.
    """
    args = parseArgs()

    folders = args.folders
    sessionFile = "session.xml"

    data = readFourierEnergyMultiple(folders)

    _, ax = plt.subplots(1, 1, figsize=(8, 6))

    for i, f in enumerate(folders):
        depth, width = extract_width_depth(f)
        time, energymodes = data[f]

        # get basefolder of f
        sessionFile = os.path.join(f, "session.xml")

        lz = extractValueXML(sessionFile, "LZ")

        # lz may be a number or some simbolic expression, like 4*depthGap, so we need to evaluate it if it's not a number
        try:
            lz = float(lz)
        except ValueError:
            # if it's not a number, we need to evaluate it
            depthGap = depth
            lz = eval(lz)
            print(
                colors.WARNING
                + f"Evaluating LZ = {lz} from {sessionFile}"
                + colors.ENDC
            )

        print(colors.OKBLUE + f"Using LZ = {lz} from {sessionFile}" + colors.ENDC)

        time = time - time[0]  # reset time to start at 0

        time, energymodes = timeFilter(time, energymodes, args.time_min, args.time_max)

        for k, energy in enumerate(energymodes):
            if k == 0 and not args.meanmode:
                continue  # skip mode 0 (mean mode)

            if k > 0:
                lz_k = lz / k
                lz_k = lz_k / depth  # nondim spatial period
            else:
                lz_k = np.inf
            if args.log:
                energy = np.log(energy)
            ax.plot(
                time,
                energy,
                label=f"Mode {k} (LZ = {lz_k:.4f}d) - d{depth}_w{width}",
                color=tableau_colors[k % len(tableau_colors)],
                linewidth=2,
                linestyle=line_styles[i % len(line_styles)],
            )
        ax.set_xlabel("t")
        ax.set_ylabel("Amplitude of Fourier Modes")

    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
