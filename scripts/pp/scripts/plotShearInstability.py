import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import os
from pp.fileManagement import readFieldsBySection, editFile
import subprocess


def addDataToPlot(dataFile: str, d: float, w: float, ax: Axes, color: str):
    xvals, yvals, data = readFieldsBySection(dataFile)

    u = data[:, :, 2]
    v = data[:, :, 3]
    u_y = data[:, :, 5]
    u_yy = data[:, :, 6]

    plotVar = u_y
        # plotVar = np.sqrt(u**2 + v**2)


    xvals_plot = xvals / w 
    # Plot the points
    for i in range(len(xvals)):
        plotVar[i] = plotVar[i] 
        if xvals[i] <= 0 or xvals[i] >= w:
            # idx_ypositive = np.where(yvals >= 0)[0][0]
            continue
        else:
            ax.plot(plotVar[i,::-1] + xvals_plot[i], yvals, color=color)


def main():
    """Main function to read data and plot points."""
    # script path
    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))
    createPointsScript = os.path.join(
        pathCurrentScript, "../../createPointsOfSectionDomainOneChkFile.sh"
    )
    lstdir = "incGapRe1000/baseflow/dns"
    _, ax = plt.subplots(figsize=(8, 6))

    # d = 3
    # w = 15
    d = 2
    w = 16
    # extension = "_stableIC"
    extension = ""

    code_name = f"d{d}_w{w}"

    # lstdir = "incNSboeingGapRe1000/directLinearSolver/blowingSuction"
    fldremote = "mesh_70_d2Udy2.fld"
    fldremote_noext, _ = os.path.splitext(fldremote)
    # we will take the baseflow.fld and the last .chk files for the modes

    # editFile(createPointsScript, replacement_line_baseflow, line_startswith)
    # subprocess.run([createPointsScript], cwd=pathCurrentScript, check=True)
    basePath = os.path.join(pathCurrentScript, "../../../src/", lstdir, code_name + extension)
    dataFile_baseflow = os.path.join(
        basePath, "data", f"points{fldremote_noext}.dat"
    )
    addDataToPlot(dataFile_baseflow, d, w, ax, "tab:blue")
    d = 4
    w = 17
    extension = "_stableIC"

    code_name = f"d{d}_w{w}"

    # lstdir = "incNSboeingGapRe1000/directLinearSolver/blowingSuction"
    fldremote = "mesh_39_d2Udy2.fld"
    fldremote_noext, _ = os.path.splitext(fldremote)
    # we will take the baseflow.fld and the last .chk files for the modes

    # editFile(createPointsScript, replacement_line_baseflow, line_startswith)
    # subprocess.run([createPointsScript], cwd=pathCurrentScript, check=True)
    basePath = os.path.join(pathCurrentScript, "../../../src/", lstdir, code_name + extension)
    dataFile_baseflow = os.path.join(
        basePath, "data", f"points{fldremote_noext}.dat"
    )

    addDataToPlot(dataFile_baseflow, d, w, ax, "tab:orange")
    ax.set_xlabel("x rescaled")
    ax.set_ylabel("y")
    ax.set_title(f"Baseflow velocity magnitude rescaled d={d}, w={w}")
    ax.grid()

    # editFile(createPointsScript, replacement_line_TSmode, line_startswith)
    # subprocess.run([createPointsScript], cwd=pathCurrentScript, check=True)
    # _, ax2 = plt.subplots(figsize=(8, 6))
    # dataFile_TSmode = os.path.join(
    #     basePath, "data", f"pointsmesh_{code_name}_avg_n{npointsYdir}.dat"
    # )
    # addDataToPlot(dataFile_TSmode, d, w, False, ax2)
    # ax2.set_xlabel("x rescaled")
    # ax2.set_ylabel("y")
    # ax2.set_title(f"TS mode velocity magnitude rescaled d={d}, w={w}")
    # ax2.grid()

    # Show the plot
    plt.show()


if __name__ == "__main__":
    main()
