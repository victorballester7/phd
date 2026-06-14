import numpy as np
import matplotlib.pyplot as plt

def main():
    # nc1 = np.loadtxt("../../data/neutral_curve_blasius_spatial.dat")
    # nc2 = np.loadtxt("../../data/neutral_curve_blasius_spatial_gaster.dat")
    nc3 = np.loadtxt("../../src/flatPlateRe1000inc/directLinearSolver/blowingSuction/wgnInsideDomainDivFree/data/neutral_curve_temporal_gaster.dat")
    nc4 = np.loadtxt("../../src/incGapRe1000/directLinearSolver/blowingSuctionCoarserMesh/d1.25_w90/data/neutral_curve_temporal_gaster.dat")
    

    _, ax = plt.subplots()

    # plt.plot(nc1[:, 0], nc1[:, 1], "o", label="spatial solver (clean Blasius)")
    # plt.plot(nc2[:, 0], nc2[:, 1], "o", label="temporal solver with gaster's tran (clean Blasius)")
    plt.plot(nc3[:, 0], nc3[:, 1], "o", label="temporal solver with gaster's trans (flat plate)")
    plt.plot(nc4[:, 0], nc4[:, 1], "o", label="temporal solver with gaster's trans (gap d=1.25, w=90)")

    ax.set_xlabel("Re_δ*")
    ax.set_ylabel("α_r")
    plt.legend()
    plt.show()





if __name__ == "__main__":
    main()
