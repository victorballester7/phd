import numpy as np
import matplotlib.pyplot as plt
from pp.colors import colors

BLASIUS_C = 1.7207876573


def main():
    cases = [
        # "d1.5_w10",
        # "d1.5_w15",
        # "d1.5_w20",
        # "d1.5_w25",
        "d1.5_w30",
        # "d1.5_w35",
        # "d1.5_w40",
        # "d1.5_w45",
        # "d1.5_w50",
        # "d0.5_w40",
        # "d0.75_w40",
        # "d1_w39",
        "d1.25_w40",
        # "d1.5_w40",
        "d1.75_w40",
        # "d1_w20",
        # "d2_w20",
        # "d2.5_w20",
        "d3_w21",
    ]
    #
    _, ax1 = plt.subplots()
    #
    # temporal = False
    # ylabel = "α_r" if temporal else "ω_r"
    #
    scatter = True

    for case in cases:
        file_path = f"/home/victor/Desktop/PhD/src/incGapRe1000/directLinearSolver/blowingSuction/{case}/data/xSections/neutral_curve_temporal_gaster.dat"
        nc = np.loadtxt(file_path)

        reynolds = nc[:, 0]
        omegas = nc[:, 1] / (reynolds * BLASIUS_C) * 1e6

        # plot with points
        if scatter:
            plt.plot(reynolds, omegas, "o", label=f"{case}")

        # plot with lines, two lines, upper and lower boundaries
        # discard the reynolds such that there is only one data point for each reynolds number, and plot the upper and lower boundaries
        else:
            for r in np.unique(nc[:, 0]):
                nc_r = nc[nc[:, 0] == r]
                if len(nc_r) != 2:
                    # remove the reynolds number from the array
                    print(
                        colors.WARNING
                        + f"Reynolds number {r} has {len(nc_r)} points for case = {case}, expected 2. Discarding this Reynolds number."
                        + colors.ENDC
                    )
                    nc = nc[nc[:, 0] != r]

            upper_boundary = nc[::2, 1]
            lower_boundary = nc[1::2, 1]

            reynolds = nc[::2, 0]
            # fix color for plotting
            col = plt.get_cmap("tab10")(cases.index(case) % 10)
            plt.plot(reynolds, upper_boundary, "-", color=col, label=f"{case}")
            plt.plot(reynolds, lower_boundary, "-", color=col)
    
    cases = [
        "d1.5_w70",
        "d1.75_w33",
        "d2_w41",
        "d3.5_w21",
        # "d1.5_w10",
        # "d1.5_w15",
        # "d1.5_w20",
        # "d1.5_w25",
        # "d1.5_w30",
        # "d1.5_w35",
        # "d1.5_w40",
        # "d1.5_w45",
        # "d1.5_w50",
        # "d0.5_w40",
        # "d0.75_w40",
        # "d1_w39",
        # "d1.25_w40",
        # "d1.5_w40",
        # "d1.75_w40",
        # "d1_w20",
        # "d2_w20",
        # "d2.5_w20",
        # "d3_w21",
    ]
    #
    #
    # temporal = False
    # ylabel = "α_r" if temporal else "ω_r"
    #
    scatter = True

    for case in cases:
        file_path = f"/home/victor/Desktop/PhD/src/incGapRe800/directLinearSolver/blowingSuction/{case}/data/xSections/neutral_curve_temporal_gaster.dat"
        nc = np.loadtxt(file_path)

        reynolds = nc[:, 0]
        omegas = nc[:, 1] / (reynolds * BLASIUS_C) * 1e6

        # plot with points
        if scatter:
            plt.plot(reynolds, omegas, "o", label=f"{case}")

        # plot with lines, two lines, upper and lower boundaries
        # discard the reynolds such that there is only one data point for each reynolds number, and plot the upper and lower boundaries
        else:
            for r in np.unique(nc[:, 0]):
                nc_r = nc[nc[:, 0] == r]
                if len(nc_r) != 2:
                    # remove the reynolds number from the array
                    print(
                        colors.WARNING
                        + f"Reynolds number {r} has {len(nc_r)} points for case = {case}, expected 2. Discarding this Reynolds number."
                        + colors.ENDC
                    )
                    nc = nc[nc[:, 0] != r]

            upper_boundary = nc[::2, 1]
            lower_boundary = nc[1::2, 1]

            reynolds = nc[::2, 0]
            # fix color for plotting
            col = plt.get_cmap("tab10")(cases.index(case) % 10)
            plt.plot(reynolds, upper_boundary, "-", color=col, label=f"{case}")
            plt.plot(reynolds, lower_boundary, "-", color=col)

    cases = [
        "d0.75_w93",
        "d1.25_w38",
        "d1.5_w33",
        "d3_w14",
        # "d1.5_w10",
        # "d1.5_w15",
        # "d1.5_w20",
        # "d1.5_w25",
        # "d1.5_w30",
        # "d1.5_w35",
        # "d1.5_w40",
        # "d1.5_w45",
        # "d1.5_w50",
        # "d0.5_w40",
        # "d0.75_w40",
        # "d1_w39",
        # "d1.25_w40",
        # "d1.5_w40",
        # "d1.75_w40",
        # "d1_w20",
        # "d2_w20",
        # "d2.5_w20",
        # "d3_w21",
    ]
    #
    #
    # temporal = False
    # ylabel = "α_r" if temporal else "ω_r"
    #
    scatter = True

    for case in cases:
        file_path = f"/home/victor/Desktop/PhD/src/incGapRe3000/directLinearSolver/blowingSuction/{case}/data/xSections/neutral_curve_temporal_gaster.dat"
        nc = np.loadtxt(file_path)

        reynolds = nc[:, 0]
        omegas = nc[:, 1] / (reynolds * BLASIUS_C) * 1e6

        # plot with points
        if scatter:
            plt.plot(reynolds, omegas, "o", label=f"{case}")

        # plot with lines, two lines, upper and lower boundaries
        # discard the reynolds such that there is only one data point for each reynolds number, and plot the upper and lower boundaries
        else:
            for r in np.unique(nc[:, 0]):
                nc_r = nc[nc[:, 0] == r]
                if len(nc_r) != 2:
                    # remove the reynolds number from the array
                    print(
                        colors.WARNING
                        + f"Reynolds number {r} has {len(nc_r)} points for case = {case}, expected 2. Discarding this Reynolds number."
                        + colors.ENDC
                    )
                    nc = nc[nc[:, 0] != r]

            upper_boundary = nc[::2, 1]
            lower_boundary = nc[1::2, 1]

            reynolds = nc[::2, 0]
            # fix color for plotting
            col = plt.get_cmap("tab10")(cases.index(case) % 10)
            plt.plot(reynolds, upper_boundary, "-", color=col, label=f"{case}")
            plt.plot(reynolds, lower_boundary, "-", color=col)
    # ax1.set_xlabel("Re_δ*")
    # ax1.set_ylabel(ylabel)
    # 
    # extra = [
    #     # "/home/victor/Desktop/PhD/data/neutralCurves/analytical_blasius_spatial.dat",
    #     # "/home/victor/Desktop/PhD/data/neutralCurves/analytical_blasius_temporal_gaster.dat.bak",
    #     "/home/victor/Desktop/PhD/data/neutralCurves/analytical_blasius_temporal_gaster.dat",
    #     "/home/victor/Desktop/PhD/src/flatPlateRe1000inc/directLinearSolver/blowingSuction/data/neutral_curve_temporal_gaster.dat",
    # ]
    # labels = ["temporal gaster", "dns temporal gaster"]
    # for i, (e, l) in enumerate(zip(extra, labels)):
    #     extra = np.loadtxt(e)
    #     reynolds = extra[:, 0]
    #     omega_r = extra[:, 1]
    #     omega_r = omega_r / (reynolds * BLASIUS_C) * 1e6
    #     plt.plot(reynolds, omega_r, "o", label=l)
    #
    # def re_to_x(re):
    #     return (
    #         ((np.asarray(re) / RE_REFERENCE) ** 2 - 1.0) * RE_REFERENCE / BLASIUS_C**2
    #     )
    #
    # def x_to_re(x):
    #     return RE_REFERENCE * np.sqrt(1.0 + BLASIUS_C**2 * np.asarray(x) / RE_REFERENCE)
    #
    # ax2 = ax1.secondary_xaxis(
    #     "top",
    #     functions=(re_to_x, x_to_re),
    # )

    # ax2.set_xlabel("x")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
