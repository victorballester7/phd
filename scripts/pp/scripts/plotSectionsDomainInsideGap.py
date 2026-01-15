import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import os
from pp.fileManagement import readFieldsBySection, editFile
import subprocess


def addDataToPlot(dataFile: str, d: float, w: float, baseflow: bool, ax: Axes):
    xvals, yvals, data = readFieldsBySection(dataFile)

    u = data[:, :, 2]
    v = data[:, :, 3]

    if baseflow:
        plotVar = u
        # plotVar = np.sqrt(u**2 + v**2)
    else:
        print("Plotting Reynolds stresses")
        uu = data[:, :, 5]
        vv = data[:, :, 7]
        plotVar = np.sqrt(np.abs(u**2 + v**2 + uu + vv))

    xvals_plot = xvals / w if baseflow else xvals / w * 4

    # Plot the points
    for i in range(len(xvals)):
        plotVar[i] = plotVar[i] / np.max(plotVar[i])
        if xvals[i] <= 0 or xvals[i] >= w:
            idx_ypositive = np.where(yvals >= 0)[0][0]
            ax.plot(
                plotVar[i, idx_ypositive:] + xvals_plot[i],
                yvals[idx_ypositive:],
                color="tab:orange",
            )
        else:
            ax.plot(plotVar[i] + xvals_plot[i], yvals, color="tab:blue")


def main():
    """Main function to read data and plot points."""
    # script path
    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))
    createPointsScript = os.path.join(
        pathCurrentScript, "../../createPointsOfSectionDomainOneChkFile.sh"
    )

    # d = 3
    # w = 15
    d = 3
    w = 26
    extension = "_stableIC"

    code_name = f"d{d}_w{w}"

    # lstdir = "incNSboeingGapRe1000/directLinearSolver/blowingSuction"
    lstdir = "incNSboeingGapRe1000/baseflow/dns"
    fldremote = f"mesh_{code_name}_25.bak0.chk"
    fldremote_noext, _ = os.path.splitext(fldremote)
    # we will take the baseflow.fld and the last .chk files for the modes

    basePath = os.path.join(pathCurrentScript, "../../../src/", lstdir, code_name + extension)

    npointsYdir = 100
    numlinesInsideGap = 30
    dxPointsInsideGap = w / numlinesInsideGap
    minY = -d
    maxY = d * 0.666
    eps = w / 5
    xloc = np.arange(-eps, w + eps + dxPointsInsideGap, dxPointsInsideGap)

    line_startswith = np.array(
        [
            "constVar=",
            "constValue=",
            "varValueMIN=",
            "varValueMAX=",
            "N=",
            "CASE=",
            "DIR=",
            "FLD_REMOTE=",
        ]
    )

    replacement_line_baseflow = np.array(
        [
            'constVar="x"',
            f"constValue=({' '.join(map(str, xloc))})",
            f"varValueMIN={minY}",
            f"varValueMAX={maxY}",
            f"N={npointsYdir}",
            f'CASE="{code_name}"',
            f'DIR="{lstdir}/{code_name}{extension}"',
            # 'FLD_REMOTE="baseflow.fld"',
            f'FLD_REMOTE="{fldremote}"',
        ]
    )

    replacement_line_TSmode = np.array(
        [
            'constVar="x"',
            f"constValue=({' '.join(map(str, xloc))})",
            f"varValueMIN={minY}",
            f"varValueMAX={maxY}",
            f"N={npointsYdir}",
            f'CASE="{code_name}"',
            f'DIR="{lstdir}/{code_name}"',
            f'FLD_REMOTE="mesh_{code_name}_avg.fld"',
        ]
    )

    _, ax = plt.subplots(figsize=(8, 6))
    # editFile(createPointsScript, replacement_line_baseflow, line_startswith)
    # subprocess.run([createPointsScript], cwd=pathCurrentScript, check=True)
    dataFile_baseflow = os.path.join(
        basePath, "data", f"points{fldremote_noext}_n{npointsYdir}.dat"
    )
    addDataToPlot(dataFile_baseflow, d, w, True, ax)
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
