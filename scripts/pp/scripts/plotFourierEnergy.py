import matplotlib.pyplot as plt
import numpy as np
from pp.fft import fftFreqs
from pp.filterData import timeFilter
from pp.inputargs import parseArgs
from pp.colors import colors
from pp.fileManagement import extract_width_depth, readFourierEnergyMultiple, extractValueXML
from pp.returnMap import getReturnPoints
import os


def main():
    """
    Main function to plot history points from one or multiple folders as a function of time or in phase space (u vs v). The points are stored in the historyPoints.dat file in each folder.
    """
    args = parseArgs()

    folders = args.folders
    sessionFile = "session.xml"

    data = readFourierEnergyMultiple(folders)

    _, ax = plt.subplots(1, 1, figsize=(8, 6))

    for f in folders:
        time, energymodes = data[f]

        # get basefolder of f
        basefolder = os.path.dirname(f)
        sessionFile = os.path.join(basefolder, "session.xml")


        lz = extractValueXML(sessionFile, "LZ")

        print(colors.OKBLUE + f"Using LZ = {lz} from {sessionFile}" + colors.ENDC)


        time = time - time[0]  # reset time to start at 0

        time, energymodes = timeFilter(time, energymodes, args.time_min, args.time_max)
        depth, width = extract_width_depth(f)

        for k, energy in enumerate(energymodes):
            if k == 0:
                continue  # skip mode 0 (mean mode)

            lz_k = lz / k
            lz_k = lz_k / depth  # nondim spatial period
            ax.plot(time, energy, label=f"Mode {k} (LZ = {lz_k:.4f}d) - d{depth}_w{width}")
        ax.set_xlabel("t")
        ax.set_ylabel("Amplitude of Fourier Modes")


    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
