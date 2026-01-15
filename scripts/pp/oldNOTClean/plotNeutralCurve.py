import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
import subprocess
from typing import Tuple
import plotNfactorOrrSommerfeld as func


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


def runOStemporalAnalysis(tomlFile: Path) -> None:
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
            "n = 175",
            "re = 1000",
            'branch = "temporal"',
            'problem = "BoundaryLayer"',
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


def getomega_rNeutralCurve(omega_r: np.ndarray, omega_i: np.ndarray) -> np.ndarray:
    # find indices where omega_i changes sign
    indices = np.where(np.diff(np.sign(omega_i)))[0]

    omegas_r_neutral = np.array([])

    for i in indices:
        idx1 = i
        idx2 = i + 1

        omega_r_neutral = (
            omega_r[idx1] * omega_i[idx2] - omega_r[idx2] * omega_i[idx1]
        ) / (omega_i[idx2] - omega_i[idx1])

        omegas_r_neutral = np.append(omegas_r_neutral, omega_r_neutral)

    return omegas_r_neutral


# def computeNfactor(a, b, alpha_i_a, alpha_i_b, nfactor_a, i):
#     trap = np.
#     # if i > 0:
#     #     trapezoidal = -(b - a) * (alpha_i_b + alpha_i_a) / 2
#     # else:
#     #     trapezoidal = -(b - a) * (alpha_i_b - alpha_i_a) / 2
#     # nfactor_b = nfactor_a + trapezoidal
#     return nfactor_b


def plot(reynolds: np.ndarray, omegas_r: np.ndarray, xpositions: np.ndarray) -> None:
    fig, ax1 = plt.subplots()

    # Plot the main data
    ax1.plot(reynolds, omegas_r, "o")
    ax1.set_xlabel("Re_δ*")
    ax1.set_ylabel("ω_r")
    ax1.grid()

    # Create the secondary x-axis (top), sharing the y-axis
    if xpositions.size > 0:
        ax2 = ax1.twiny()
        
        # Match x-axis limits of the bottom axis
        ax2.set_xlim(ax1.get_xlim())

        # Choose ticks (can be customized: e.g. np.linspace, or just use all values)
        tick_locs = reynolds
        tick_labels = [f"{x}" for x in xpositions]

        ax2.set_xticks(tick_locs)
        ax2.set_xticklabels(tick_labels)
        ax2.set_xlabel("x")

    plt.title("Neutral curve ω_i = 0")
    plt.tight_layout()
    plt.show()


def runAnalyticalBlasius(
    re_vary_bacause_of_deltaStar: bool,
    tomlFile: Path,
    filename_ev: Path,
    alpha_r_min: float,
    alpha_r_max: float,
    alpha_r_num: int,
) -> None:
    omegas_r = np.array([])

    reynolds = np.array([])
    xpositions = np.array([])

    if re_vary_bacause_of_deltaStar:
        # don't go beyond -150, because the tip of the neutral curve is around there
        xpositions = np.arange(-150, 1001, 50)

        C = 1.7207876573
        alpha_r_min = alpha_r_min / 2
        for x in xpositions:
            delta_star_x = np.sqrt(1 + x * C**2 / 1000)
            re = 1000 * delta_star_x

            reynolds = np.append(reynolds, re)

            normalization_factor = delta_star_x

            print(f"Running runs for x = {x} and Re = {re}")

            replacement_line = np.array(
                [
                    f"re = {re * normalization_factor}",
                    'problem = "BoundaryLayer"',
                    f'filenameUprofile = "/home/victor/Desktop/PhD/src/flatSurfaceRe1000IncNS/linearSolver_blowingSuction/data/points_x{x}_n800.dat"',
                    f"vars_r = {{min = {alpha_r_min * normalization_factor}, max = {alpha_r_max * normalization_factor}, num = {alpha_r_num}}}",
                ]
            )
            line_startswith = np.array(
                [
                    "re = ",
                    "problem = ",
                    "filenameUprofile = ",
                    "vars_r = ",
                ]
            )

            # Read and modify the toml file
            editFile(tomlFile, replacement_line, line_startswith)

            runOStemporalAnalysis(tomlFile)
            omega_r, omega_i, _ = func.readData(filename_ev)

            omega_r = omega_r / normalization_factor
            omega_i = omega_i / normalization_factor

            o_r = getomega_rNeutralCurve(omega_r, omega_i)

            omegas_r = np.append(omegas_r, o_r)
    else:
        re_min = 520
        re_max = 2000
        re_num = 50
        reynolds = np.linspace(re_min, re_max, re_num)
        for re in reynolds:
            print(f"Running runs for Re = {re}")

            replacement_line = np.array([f"re = {re}"])
            line_startswith = np.array(["re = "])
            # Read and modify the toml file
            editFile(tomlFile, replacement_line, line_startswith)

            runOStemporalAnalysis(tomlFile)
            omega_r, omega_i, _ = func.readData(filename_ev)

            o_r = getomega_rNeutralCurve(omega_r, omega_i)

            omegas_r = np.append(omegas_r, o_r)

    reynolds = np.repeat(reynolds, 2)  # there are two values of omega_r for each Re
    xpositions = np.repeat(xpositions, 2)  # there are two values of omega_r for each Re
    print(reynolds)
    print(omegas_r)
    omegas_r = omegas_r.flatten()
    plot(reynolds, omegas_r,xpositions=xpositions)


def runSectionsBaseflowDNS(tomlFile: Path, filename_ev: Path) -> None:
    xpositions = np.arange(100, 1001, 50)
    # add -70 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 70 80 90 
    # xpositions = np.insert(xpositions, 0, np.arange(-70, 0, 10))
    # xpositions = np.insert(xpositions, 0, np.arange(10, 100, 10))

    omegas_r = np.array([])
    omega_r_xpositions = []
    C = 1.7207876573
    n_interp_dns = 600  # number of points in the DNS interpolation

    for x in xpositions:
        re = 1000 * np.sqrt(1 + x * C**2 / 1000)
        print(f"Running runs for x = {x} and Re = {re}")

        replacement_line = np.array(
            [
                f"re = {re}",
                'problem = "Custom"',
                f'filenameUprofile = "/home/victor/Desktop/PhD/src/incNSboeingGapRe1000/baseflow/dns/d2_w36/data/points90_x{x}_n{n_interp_dns}.dat"'
            ]
        )
        line_startswith = np.array(["re = ", "problem = ", "filenameUprofile = "])

        # Read and modify the toml file
        editFile(tomlFile, replacement_line, line_startswith)

        runOStemporalAnalysis(tomlFile)
        omega_r, omega_i, _ = func.readData(filename_ev)

        o_r = getomega_rNeutralCurve(omega_r, omega_i)

        omegas_r = np.append(omegas_r, o_r)
        omega_r_xpositions.append([int(x), o_r.tolist()])

    xpositions = np.repeat(xpositions, 2)  # there are two values of omega_r for each Re
    print(omega_r_xpositions)
    omegas_r = omegas_r.flatten()

    reynolds = 1000 * np.sqrt(1 + xpositions * C**2 / 1000)

    plot(reynolds, omegas_r,xpositions=xpositions)


def main():
    # PARAMETERS TO CHANGE
    alpha_r_min = 0.01
    alpha_r_max = 0.5
    alpha_r_num = 50

    analyticalBlasius = False

    # if analyticalBlasius = True
    re_vary_bacause_of_deltaStar = True  # in this case (same as in the DNS simulations), we need to propperly rescale Re alpha and the outpu omega_r

    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))
    # xpositions = getDatafromDNS(n_interp_dns, pathCurrentScript)

    tomlFile = "../../../../orrSommerfeld/config/input.toml"
    tomlFile = Path(os.path.join(pathCurrentScript, tomlFile))
    filename_ev = "../../../../orrSommerfeld/data/eigenvalues.dat"
    filename_ev = Path(os.path.join(pathCurrentScript, filename_ev))

    # Setup the toml file
    setupToml(tomlFile, alpha_r_min, alpha_r_max, alpha_r_num)

    if analyticalBlasius:
        runAnalyticalBlasius(
            re_vary_bacause_of_deltaStar,
            tomlFile,
            filename_ev,
            alpha_r_min,
            alpha_r_max,
            alpha_r_num,
        )
    else:
        runSectionsBaseflowDNS(tomlFile, filename_ev)


if __name__ == "__main__":
    main()
