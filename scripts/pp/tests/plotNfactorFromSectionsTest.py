import numpy as np
import matplotlib.pyplot as plt
import os
import re
from matplotlib.axes import Axes


def read_points_from_file(filename):
    """Reads points from a file and returns them as a numpy array."""
    return np.loadtxt(filename, skiprows=3)


def addNxplot(dataFile: np.ndarray, avgFiles: bool, ax: Axes):
    """Main function to read data and plot points."""

    points_tmp = []

    rms = []
    var = np.array([])

    data = []

    xvals = np.array([])
    yvals = np.array([])

    if avgFiles:
        for dF in dataFile:
            data_dF = read_points_from_file(dF)
            data_tmp = []
            if len(xvals) == 0:
                xvals = np.unique(data_dF[:, 0])
            if len(yvals) == 0:
                yvals = np.unique(data_dF[:, 1])
                yvals.sort()
            for x in xvals:
                tmp = data_dF[data_dF[:, 0] == x]
                tmp = tmp[:, [2, 3]]
                data_tmp.append(tmp)  # we only store u and v
            data.append(data_tmp)

        data = np.array(data)

        rms = np.sqrt(np.mean(data[:, :, :, 0] ** 2 + data[:, :, :, 1] ** 2, axis=0))
        ymin = 0.2
        yindx = np.argmin(np.abs(yvals - ymin))
        Loo = np.max(np.abs(rms[:, yindx:]), axis=1)
        ymax = 150
        yindx = np.argmin(np.abs(yvals - ymax))
        L2 = np.sqrt(np.trapezoid(rms[:, :yindx] ** 2, yvals[:yindx]))
        xvals = xvals[1:]

        Aoo, A2 = Loo[1:], L2[1:]  # skip the first point
        Nxoo = np.log(Aoo / Aoo[0])
        Nx2 = np.log(A2 / A2[0])
        ax.plot(xvals, Nxoo, markersize=2, label=f"{dataFile[0]} (Loo)", linestyle="--")
        ax.plot(xvals, Nx2, markersize=2, label=f"{dataFile[0]} (L2)", linestyle="-")
    else:
        for dF in dataFile:
            data_dF = read_points_from_file(dF)
            data_tmp = []
            if len(xvals) == 0:
                xvals = np.unique(data_dF[:, 0])
            if len(yvals) == 0:
                yvals = np.unique(data_dF[:, 1])
                yvals.sort()
            for x in xvals:
                tmp = data_dF[data_dF[:, 0] == x]
                data_tmp.append(tmp) # we only store u and v
            data.append(data_tmp)

        data = np.array(data)

        u = data[:, :, :, 2]
        v = data[:, :, :, 3]
        uu = data[:, :, :, 5]
        vv = data[:, :, :, 7]
        rms = np.sqrt(np.abs(u**2 + v**2 + uu + vv))[0]
        # rms = np.sqrt(np.abs(data[:, :, :, -1]))[0]
        Loo = np.max(np.abs(rms), axis=1)
        ymax = 150
        yindx = np.argmin(np.abs(yvals - ymax))
        L2 = np.sqrt(np.trapezoid(rms[:, :yindx] ** 2, yvals[:yindx]))
        LooArgmin = np.argmin(Loo)
        L2Argmin = np.argmin(L2)
        Aoo, A2 = Loo[LooArgmin:], L2[L2Argmin:]  # skip to the minimum point
        xvals = xvals[LooArgmin:]
        Nxoo = np.log(Aoo / Aoo[0])
        Nx2 = np.log(A2 / A2[0])
        ax.plot(xvals, Nxoo, markersize=2, label=f"{dataFile[0]} (Loo)", linestyle="--")
        ax.plot(xvals, Nx2, markersize=2, label=f"{dataFile[0]} (L2)", linestyle="-")
    


def get_chkfiles_numbers(basePath: str, n: int) -> np.ndarray:
    """Get all the chk numbers of the files of the form points{chk}_all_n{n}.dat inside basePath."""
    pattern = re.compile(rf"points(\d+)_all_n{n}\.dat")
    chkfiles = []
    for filename in os.listdir(basePath):
        match = pattern.match(filename)
        if match:
            chkfiles.append(int(match.group(1)))
    chkfiles.sort()
    return np.array(chkfiles)


def getDataFileList(basePath: str, n: int) -> np.ndarray:
    """Get all the data files of the form points{chk}_all_n{n}.dat inside basePath."""
    chkfiles = get_chkfiles_numbers(basePath, n)
    dataFile = [f"points{chk}_all_n{n}.dat" for chk in chkfiles]
    dataFile = [os.path.join(basePath, f) for f in dataFile]
    return np.array(dataFile)


def main():
    # script path
    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))

    basePath = "../src/flatSurfaceRe1000IncNS/directLinearSolver/blowingSuction" 
    basePath = os.path.join(pathCurrentScript, basePath)

    fig, ax = plt.subplots(figsize=(8, 6))

    folders1 = [
        "wgnInsideDomainDivFree",
    ]
    n = 600
    for f in folders1:
        basePath_dw = os.path.join(basePath, f, "data")
        dataFile = np.array([os.path.join(basePath_dw, f"pointsavg_all_n{n}.dat")])
        addNxplot(dataFile, False, ax)

    basePath = "../src/incNSboeingGapRe1000/directLinearSolver/blowingSuction/d4_w15"
    basePath = os.path.join(pathCurrentScript, basePath)

    folders2 = [
        # "wgnInsideDomainNektarAveraged",
        # "wgnInsideDomainOnlyVNektarAveraged",
        # "wgnInsideDomainSameMagnitudeNektarAveraged",
        # "wgnInsideDomainOnlyVy0.1NektarAveraged",
        # "wgnInsideDomainOnlyVy1NektarAveraged",
        "wgnInsideDomainDivFreeNektarAveraged",
    ]

    n = 600
    for f in folders2:
        basePath_dw = os.path.join(basePath, f, "data")
        dataFile = np.array([os.path.join(basePath_dw, f"pointsavg_all_n{n}.dat")])
        addNxplot(dataFile, False, ax)

    ax.legend()
    ax.grid()

    plt.show()


if __name__ == "__main__":
    main()
