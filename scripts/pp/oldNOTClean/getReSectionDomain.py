import numpy as np
from src.blasiusIncNS.physicalQuantities import computeDelta
import os
import sys
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

DELTASTAR = 1.7207876575203542
DELTASTARtoDELTA = 4.910287765963463 / DELTASTAR


def readData(fileName: str):
    data = np.loadtxt(fileName, skiprows=3)
    y = data[:, 1]
    u = data[:, 2]

    return y, u

def computeDeltastar_dns(filename: str, x_values: np.ndarray):
    deltaStar_dns = []
    # get directory of the current file
    pwd = os.path.dirname(os.path.abspath(__file__))

    for x in x_values:
        filename_x = filename + f"{x}.dat"
        path = os.path.join(pwd, filename_x)

        y, u = readData(path)

        delta = computeDelta(y, u)

        deltaStar = delta / DELTASTARtoDELTA
        nu = 1.0 / 1000
        U = 1.0
        re_deltaStar = deltaStar * U / nu
        deltaStar_dns.append(deltaStar)
        print("x = ", x)
        print("   delta = ", delta)
        print("   Re_deltaStar = ", re_deltaStar)

    return deltaStar_dns

def main():
    # get directory of the current file
    pwd = os.path.dirname(os.path.abspath(__file__))
    filename = "../src/flatSurfaceRe1000IncNS/dns/changedformula_newUpperBC/data"
    filename += "/points_x"
    x_values = np.array(
        [
            -50,
            -30,
            -10,
            0,
            50,
            100,
            150,
            200,
            250,
            300,
            350,
            400,
            450,
            500,
            550,
            600,
            650,
            700,
            750,
            800,
            850,
            900,
            950,
            1000,
        ]
    )
    deltaStar_dns = []
    deltastar_theo = np.sqrt(1 + x_values * DELTASTAR**2 / 1000)

    deltaStar_dns = computeDeltastar_dns(filename, x_values)

    # filename = "../src/flatSurfaceRe1000IncNS/dns/data"
    # filename += "/points_x"

    # deltaStar_dns_old =computeDeltastar_dns(filename, x_values)


    # plot the results
    plt.figure()
    # x_values = np.log10(x_values)
    # deltaStar_dns = np.log10(deltaStar_dns)
    # deltastar_theo = np.log10(deltastar_theo)
    # plt.plot(x_values, deltaStar_dns, label="deltaStar_dns")
    # plt.plot(x_values, deltastar_theo, label="deltaStar_theo")
    plt.plot(x_values, np.abs(deltaStar_dns - deltastar_theo)/np.abs(deltastar_theo), label="|deltaStar_dns - deltaStar_theo|/deltaStar_theo")

    plt.xlabel("x")
    plt.ylabel("deltaStar")
    plt.legend()
    plt.grid()

    plt.show()


if __name__ == "__main__":
    main()
