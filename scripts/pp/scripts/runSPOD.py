import numpy as np
import matplotlib.pyplot as plt
import os
from pp.fileManagement import extract_depth_width, readDataHistoryPoints
from pp.inputargs import parseArgs
from pp.filterData import timeFilter, spatialFilter
from pp.spod import plot_confidence_bounds, spod, plot_spectrum
import h5py

def getFreqMode1(L, f):
    L_mode1 = L[:, 0]  # first mode
    idx_peak = np.argmax(L_mode1)

    # iterpolate with parabola to find more accurate peak frequency
    if idx_peak == 0 or idx_peak == len(L_mode1) - 1:
        f_peak = f[idx_peak]
    else:
        y0 = L_mode1[idx_peak - 1]
        y1 = L_mode1[idx_peak]
        y2 = L_mode1[idx_peak + 1]
        x0 = f[idx_peak - 1]
        x1 = f[idx_peak]
        x2 = f[idx_peak + 1]

        # Parabolic interpolation formula
        f_peak = x1 - 0.5 * ((y2 - y0) / (y2 - 2 * y1 + y0)) * (x2 - x0)
    print(f"Peak frequency (ω) of mode 1: {f_peak * 2*np.pi:.6f}")


def main():
    # get path working directory
    args = parseArgs()

    folder = args.folders[0]  # we do one by one
    time_min = args.time_min
    time_max = args.time_max

    d, w = extract_depth_width(folder)

    xmin = 0
    xmax = 3 * w
    ymin = -d
    ymax = 2

    points, time, fields = readDataHistoryPoints(folder)
    points, fields = spatialFilter(points, fields, xmin, xmax, ymin, ymax)
    time, fields = timeFilter(time, fields, time_min, time_max)

    # remove pressure field (last variable)
    fields = fields[:, :, :-1]
    print(f"Data shape after filtering: {fields.shape}")

    # convert 22 * 13400 * 2 array to (13400, 22*2) matrix (nrows = time steps, ncols = points * vars)
    fields = np.reshape(
        np.swapaxes(fields, 0, 1), (fields.shape[1], points.shape[0] * fields.shape[2])
    )
    # make time mean zero
    # fields = fields - np.mean(fields, axis=0)
    dt = time[1] - time[0]

    if not os.path.isdir(folder):
        folder = os.path.dirname(folder)

    # apparently the frequencies do not depend on the spatial weight function
    spod(
        fields,
        dt,
        folder,
        weight="default",
        nOvlp="default",
        # nDFT=2000,
        nDFT="default",
        window="default",
        method="fast",
    )

    SPOD_LPf = h5py.File(
        os.path.join(folder, "SPOD_LPf.h5"), "r"
    )  # load data from h5 format
    L = SPOD_LPf["L"][:, :]  # modal energy E(f, M)
    P = SPOD_LPf["P"][:, :, :]  # mode shape
    f = SPOD_LPf["f"][:]  # frequency
    SPOD_LPf.close()

    getFreqMode1(L, f)

    plot_spectrum(f,L,hl_idx=5,loglog=False)
    plt.show()


if __name__ == "__main__":
    main()
