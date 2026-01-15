import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import os
from pp.fileManagement import readFieldsBySection
from pp.filterData import getRMS
from pp.codenamesNfactor import code_names
from pp.parameters import Parameters


def addDataToPlot(dataFile: str, ax: Axes):
    xvals, yvals, data = readFieldsBySection(dataFile)

    rms = getRMS(data)

    # Plot the points
    for i in range(len(xvals)):
        ax.plot(rms[i], yvals, "-", markersize=2, label=f"x = {xvals[i]}")


def main():
    """Main function to read data and plot points."""
    # script path
    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))

    # dataFile = "../../../src/flatSurfaceRe1000Ma0.6ComNS/dns_shortDomain/data/pointsmesh_flat_180_n100.dat"
    dataFile = (
        "../../../src/flatSurfaceRe1000IncNS/dns/data/pointsmesh_flat_45_n100.dat"
    )
    incNS = True if "inc" in dataFile.lower() else False

    dataFile = os.path.join(pathCurrentScript, dataFile)

    fig, ax = plt.subplots(figsize=(8, 6))

    xvals, yvals, data = readFieldsBySection(dataFile)

    u = data[:, :, 2] if incNS else data[:, :, 6]
    deltas = []
    for i in range(len(xvals)):
        d = Parameters.computeDelta(yvals, u[i, :])
        deltas.append(d)
    deltas = np.array(deltas)

    if incNS:
        C_dstar = 1.7012012016016654
        C_d = 4.910088933702247
    else:
        # for Ma = 0.6
        C_dstar = 1.7610146771234572
        C_d = 4.90721172020025

    re_dstar = 1000.0
    delta_theo = C_d / C_dstar * np.sqrt(1 + xvals * C_dstar**2 / re_dstar)

    ax.plot(xvals, deltas, "-o", markersize=4, label="dns")
    ax.plot(xvals, delta_theo, "-s", markersize=4, label="theoretical")

    ax2 = ax.twinx()
    ax2.plot(
        xvals,
        np.abs(deltas - delta_theo) / np.abs(delta_theo),
        "-x",
        markersize=4,
        label="relative error",
    )

    # Set labels and title
    ax.set_xlabel("x")
    ax.set_ylabel("δ_99")
    ax2.set_ylabel("Relative error")
    ax.set_title(
        "Boundary Layer Thickness Comparison" + (" (incNS)" if incNS else " (comNS)")
    )

    # Add grid and legend
    ax.grid()
    ax.legend()

    # Show the plot
    plt.show()


if __name__ == "__main__":
    main()
