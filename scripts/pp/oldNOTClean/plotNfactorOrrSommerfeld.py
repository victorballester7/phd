import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
import subprocess
from typing import Tuple


def readData(filename: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(filename, skiprows=1)

    omega_r = data[:, 0]
    omega_i = data[:, 1]
    alpha_r = data[:, 2]

    return omega_r, omega_i, alpha_r


def getMaxTemporalAnalysis(
    omega_r: np.ndarray, omega_i: np.ndarray, alpha_r: np.ndarray
) -> int:
    # Find the index of the maximum value of omega_i (most unstable mode)
    max_index = int(np.argmax(omega_i))

    return max_index


def domega_rdalpha_r(omega_r: np.ndarray, alpha_r: np.ndarray, max_index: int) -> float:
    # Calculate the derivative of omega_r with respect to alpha_r using non-uniform grid (just in case)

    h0 = alpha_r[max_index] - alpha_r[max_index - 1]
    h1 = alpha_r[max_index + 1] - alpha_r[max_index]

    diff = (
        omega_r[max_index + 1] * h0 / (h1 * (h0 + h1))
        + omega_r[max_index] * (h1 - h0) / (h0 * h1)
        - omega_r[max_index - 1] * h1 / (h0 * (h0 + h1))
    )
    return diff


def gasterTransform(
    omega_r: np.ndarray, omega_i: np.ndarray, alpha_r: np.ndarray
) -> float:
    max_index = getMaxTemporalAnalysis(omega_r, omega_i, alpha_r)
    diff = domega_rdalpha_r(omega_r, alpha_r, max_index)

    omega_i_temp = omega_i[max_index]
    print(f"α_r_temp = {alpha_r[max_index]}")
    print(f"ω_r_temp = {omega_r[max_index]}")
    print(f"ω_i_temp = {omega_i_temp}")
    print(f"d(ω_r)/d(α_r) = {diff}")

    alpha_i = -omega_i_temp / diff

    print(f"α_i_spat = {alpha_i}")

    return alpha_i


def editFile(pathFile: Path, replacement_line, line_startswith) -> None:
    # Read the file and replace the specified line
    with open(pathFile, "r") as file:
        lines = file.readlines()

    with open(pathFile, "w") as file:
        for line in lines:
            linechanged = False
            # loop in the i,j in (replacement_line, line_startswith)
            for rep_line, line_start in zip(replacement_line, line_startswith):
                if line.strip().startswith(line_start):
                    # take the last part of the line (starting with # ) and add it to rep_line
                    rep_line = rep_line + " #" + line.split("#")[-1]
                    file.write(rep_line)
                    linechanged = True
                    break
            if not linechanged:
                file.write(line)


def getDatafromDNS(n_interp_dns, pathCurrentScript: str) -> np.ndarray:
    # Define the x positions
    xpositions = np.arange(50, 1001, 50)
    xpositions = np.insert(xpositions, 0, 20)

    # Create the string to replace the X= line
    values_str = " ".join(str(x) for x in xpositions)
    replacement_line = np.array([f"X=({values_str})", f"N={n_interp_dns}"])
    line_startswith = np.array(["X=(", "N="])

    # Path to the script
    filenameScript = "../scripts/createPointsOfSectionDomain.sh"
    filenameScript = Path(os.path.join(pathCurrentScript, filenameScript))

    # Read and modify the script
    editFile(filenameScript, replacement_line, line_startswith)

    # Make the script executable and run it
    subprocess.run([filenameScript], check=True)

    return xpositions


def runOStemporalAnalysis(uprof_filename: str, tomlFile: Path) -> None:
    replacement_line = np.array([uprof_filename])
    line_startswith = np.array(["filenameUprofile = "])

    # Read and modify the toml file
    editFile(tomlFile, replacement_line, line_startswith)

    dirMakefile = Path(os.path.join(os.path.dirname(tomlFile), "../"))

    # run 'make run' in the dirMakefile directory
    subprocess.run(["make", "run"], cwd=dirMakefile, check=True)

    return


def setupToml(
    filenameToml: Path, alpha_r_min: float, alpha_r_max: float, alpha_r_num: int
) -> None:
    line_startswith = np.array(
        [
            "n = ",
            "re = ",
            "branch = ",
            "problem = ",
            "doPlot = ",
            "use_c = ",
            "multipleRun = ",
            "plotUprofile = ",
            "vars_r = ",
            "vars_i = ",
        ]
    )
    replacement_line = np.array(
        [
            "n = 200",
            "re = 1000",
            'branch = "temporal"',
            'problem = "Custom"',
            "doPlot = false",
            "use_c = false",
            "multipleRun = true",
            "plotUprofile = false",
            f"vars_r = {{min = {alpha_r_min}, max = {alpha_r_max}, num = {alpha_r_num}}}",
            "vars_i = {min = 0.0, max = 0.0, num = 1}",
        ]
    )
    editFile(filenameToml, replacement_line, line_startswith)

    return


# def computeNfactor(a, b, alpha_i_a, alpha_i_b, nfactor_a, i):
#     trap = np.
#     # if i > 0:
#     #     trapezoidal = -(b - a) * (alpha_i_b + alpha_i_a) / 2
#     # else:
#     #     trapezoidal = -(b - a) * (alpha_i_b - alpha_i_a) / 2
#     # nfactor_b = nfactor_a + trapezoidal
#     return nfactor_b


def plot(xpositions, alphas_i, nfactor):
    fig, ax1 = plt.subplots()

    # First axis for nfactor
    ax1.plot(xpositions, nfactor, marker="o", color="tab:blue", label="n(x)")
    ax1.set_xlabel("x")
    ax1.set_ylabel("n(x)", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")
    ax1.grid()

    # Second axis for alphas_i
    ax2 = ax1.twinx()
    ax2.plot(xpositions, alphas_i, marker="o", color="tab:orange", label="-α_i(x)")
    ax2.set_ylabel("-α_i(x)", color="tab:orange")
    ax2.tick_params(axis="y", labelcolor="tab:orange")

    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="best")

    plt.title("n(x) and -α_i(x) vs x")
    plt.show()


def main():
    # PARAMETERS TO CHANGE
    widthGap = 15
    n_interp_dns = 800
    alpha_r_min = 0.05
    alpha_r_max = 0.4
    alpha_r_num = 50
    do_analysis = True

    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))
    # xpositions = getDatafromDNS(n_interp_dns, pathCurrentScript)
    xpositions = np.arange(0, 1001, 50)
    xpositions = np.insert(xpositions, 1, 20)
    xpositions = np.insert(xpositions, 0, -25)
    xpositions = np.insert(xpositions, 0, -50)
    xpositions = np.insert(xpositions, 0, -70)

    tomlFile = "../../orrSommerfeld/config/input.toml"
    tomlFile = Path(os.path.join(pathCurrentScript, tomlFile))
    filename_ev = "../../orrSommerfeld/data/eigenvalues.dat"
    filename_ev = Path(os.path.join(pathCurrentScript, filename_ev))

    if do_analysis:
        alphas_i = np.array([])
        xs = np.array([])
        nfactor = np.array([0.0])

        # Setup the toml file
        setupToml(tomlFile, alpha_r_min, alpha_r_max, alpha_r_num)

        for x in xpositions:
            # get index of x in xpositions
            xs = np.append(xs, x)

            print(f"Computing N-factor for x = {x}")
            # uprof_filename = f'filenameUprofile = "/home/victor/Desktop/PhD/src/boeingGapRe1000IncNS/baseflow/dns/d4_w{widthGap}/data/points_x{x}_n{n_interp_dns}.dat"'
            # uprof_filename = f"/home/victor/Desktop/PhD/src/incNSboeingGapRe1000/directLinearSolver/blowingSuction/d1.5_w45/wgn/data/points_x{x}_n{n_interp_dns}.dat"
            uprof_filename = f"/home/victor/Desktop/PhD/src/flatSurfaceRe1000IncNS/directLinearSolver/blowingSuction/data/points_x{x}_n{n_interp_dns}.dat"
            uprof_filename = f'filenameuprofile = "{uprof_filename}"'
            runOStemporalAnalysis(uprof_filename, tomlFile)
            omega_r, omega_i, alpha_r = readData(filename_ev)

            alpha_i_x = gasterTransform(omega_r, omega_i, alpha_r)

            alphas_i = np.append(alphas_i, alpha_i_x)
            if len(xs) == 1:
                continue

            n_x = np.trapezoid(-alphas_i, xs)
            nfactor = np.append(nfactor, n_x)

        plot(xpositions, -alphas_i, nfactor)
    else:
        # just compute the gaster transform for the already created data
        omega_r_temp, omega_i_temp, alpha_r_temp = readData(filename_ev)
        gasterTransform(omega_r_temp, omega_i_temp, alpha_r_temp)


if __name__ == "__main__":
    main()
