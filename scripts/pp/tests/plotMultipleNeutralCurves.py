import numpy as np
import matplotlib.pyplot as plt
from pp.colors import colors

BLASIUS_C = 1.7207876573


def main():
    cases = [
        "d0.5_w90",
        "d3_w16",
    ]
    rey = ["800", "1000", "3000"]
    _, ax1 = plt.subplots()

    for r in rey:
        for case in cases:
            file_path = f"/home/victor/Desktop/PhD/src/incGapRe{r}/directLinearSolver/blowingSuction/{case}/data/xSections/neutral_curve_temporal_gaster.dat"
            try:
                nc = np.loadtxt(file_path)
            except Exception as e:
                print(
                    colors.WARNING
                    + f"Could not load file {file_path}. Skipping case {case} for Re = {r}. Error: {e}"
                    + colors.ENDC
                )
                continue

            reynolds = nc[:, 0]
            omegas = nc[:, 1] / (reynolds * BLASIUS_C) * 1e6
            x = nc[:, 3]
            reynolds *= BLASIUS_C

            plt.plot(x, omegas, "o", label=f"Re = {r}, {case}")

    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
