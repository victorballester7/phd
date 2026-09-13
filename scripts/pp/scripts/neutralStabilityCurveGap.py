
from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import os
from pp.fileManagement import readFieldsBySection
from pp.filterData import getRMS
from pp.codenamesNfactor import code_names
from pp.DeltaN_computation import computeAmplitude

def getGrowthRate(x: np.ndarray, A: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    # we assume that the values of x are thought to compute the derivative every two points, like -71 -69 -61 -59 -51 -49 -41 -39 -31 -29 -21...
    
    x_plus_eps = x[1::2]
    x_minus_eps = x[::2]
    A_plus_eps = A[1::2]
    A_minus_eps = A[::2]

    growth_rate = np.log(A_plus_eps / A_minus_eps) / (x_plus_eps - x_minus_eps)
    x_mid = (x_plus_eps + x_minus_eps) / 2

    return x_mid, growth_rate


def main():
    """Main function to read data and plot points."""
    # script path
    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))

    basePath = (
        # "../../../src/incGapRe1000/directLinearSolver/omegaBlowSuct/"
        "../../../src/flatPlateRe1000inc/directLinearSolver/omegaBlowSuct/"
    )
    basePath = os.path.join(pathCurrentScript, basePath)

    code_names2 = [
        # "d2_w24/omega0.04",
        # "d2_w24/omega0.08",
        # "d2_w24/omega0.12",
        # "d2_w24/omega0.16",
        "omega0.02",
        "omega0.04",
        "omega0.1",
        "omega0.14",
        "omega0.19",
    ]
    n = 600
    chkfile = "avg"

    _, ax = plt.subplots(figsize=(8, 6))


    for dw in code_names2:
        dataFile_dw = os.path.join(basePath, dw, "data", f"points{chkfile}_n{n}.dat")
        x, A = computeAmplitude(dataFile_dw, doLoo=False, field="rms")
        x, alpha = getGrowthRate(x, A)
        ax.plot(x, alpha, "-", label=f"{dw}")


    dataFile_dw = os.path.join(basePath, "omega0.04", "data", f"points{chkfile}_n{n}old.dat")
    x, A = computeAmplitude(dataFile_dw, doLoo=False, field="rms")
    x, alpha = getGrowthRate(x, A)
    ax.plot(x, alpha, "-", label="omega=0.04 old")

    # dataFile = os.path.join(basePath, f"d2_w24/omega0.08/data/pointsavg_n{n}_y-2.dat")
    # x, A = computeAmplitude(dataFile, doLoo=False)
    # x, alpha = getGrowthRate(x, A)
    # ax.plot(x, alpha, "-", label="omega=0.08 ymin-2")
  

    # Set labels and title
    ax.set_xlabel("x")
    ax.set_ylabel("alpha")

    # Add grid and legend
    ax.grid()
    ax.legend()

    # Show the plot
    plt.show()


if __name__ == "__main__":
    main()
