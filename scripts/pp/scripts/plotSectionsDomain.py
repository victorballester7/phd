import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import os
from pp.fileManagement import readFieldsBySection
from pp.filterData import getRMS
from pp.codenamesNfactor import code_names

def addDataToPlot(dataFile: str, ax: Axes):
    xvals, yvals, data = readFieldsBySection(dataFile)

    rms = getRMS(data)


   #  # this is temporary, comment to plot reynold stresses
   #  u = data[:, :, 2]
   #  v = data[:, :, 7]
   #  rms = u
   # 
   #  etaMax = 9
   #  yMax = 4.26429554816 

   #  def eta(y):
   #      if y < yMax:
   #          return 2.120192718897813 * (y + a1_y2eta * np.tanh(b1_y2eta * y) + a2_y2eta * np.tanh(b2_y2eta * y) ** 2)
   #      else:
   #          return etaMax

   #  def rhofit(y, a_coeffs, b_coeffs):
   #      e = eta(y)
   #      num = sum(a * e**(i+1) for i, a in enumerate(a_coeffs))
   #      den = 1 + sum(b * e**(i+1) for i, b in enumerate(b_coeffs))
   #      return 1 + num / den

   #  def ufit(y, a_coeffs, b_coeffs):
   #      e = eta(y)
   #      num = sum(a * e**(i+1) for i, a in enumerate(a_coeffs))
   #      den = 1 + sum(b * e**(i+1) for i, b in enumerate(b_coeffs))
   #      return num / den

   #  rho = np.array([rhofit(y, a_rho, b_rho) for y in yvals])
   #  rhou = np.array([ufit(y, a_rhou, b_rhou) for y in yvals])
   #  
   #  ax.plot(rho, yvals, "k--", label="u profile")
   #  # ax.plot(eta(yvals), yvals, "r--", label="eta profile")
   #  
   #  # yvals = yvals / yMax * etaMax

   #  # indx_eta = np.argmin(np.abs(yvals - 12))
   #  indx_eta = -1
   #  yvals = yvals[:indx_eta]
   #  rms = rms[:, :indx_eta]
   #  ####################3

    # Plot the points
    for i in range(len(xvals)):
        ax.plot(rms[i], yvals, "-", markersize=2, label=f"x = {xvals[i]}")


def main():
    """Main function to read data and plot points."""
    # script path
    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))

    basePath = (
        # "../../../src/incNSboeingGapRe1000/directLinearSolver/blowingSuction/"
        "../../../src/flatSurfaceRe1000Ma0.6ComNS/dns/"
    )
    basePath = os.path.join(pathCurrentScript, basePath)

    code_names2 = [
        "",
        # "d3.75_w10",
    ]
    n = 800
    chkfile = "13"
    # chkfile = "avg"

    fig, ax = plt.subplots(figsize=(8, 6))


    for dw in code_names2:
        dataFile_dw = os.path.join(basePath, dw, "data", f"points{chkfile}_all_n{n}.dat")
        addDataToPlot(dataFile_dw, ax)

    # Set labels and title
    ax.set_xlabel("x")
    ax.set_ylabel("sqrt(u^2 + v^2)")

    # Add grid and legend
    ax.grid()
    ax.legend()

    # Show the plot
    plt.show()


if __name__ == "__main__":
    main()
