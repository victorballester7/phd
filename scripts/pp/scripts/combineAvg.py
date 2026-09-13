import os
from pp.ssh.fieldconvert import combineAvg


if __name__ == "__main__":
    # case = "d0.25_w40"
    # dir = f"/home/victor/Desktop/PhD/src/incGapRe1000/directLinearSolver/blowingSuction/{case}/"

    # dir = "/home/victor/Desktop/PhD/src/flatPlateRe800inc/directLinearSolver/blowingSuction/"

    dir = "/home/victor/Desktop/PhD/src/bfsRe1000inc/directLinearSolver/blowingSuction/d1.5/"
    combineAvg(dir, chkfile_ending_pattern="_stress.fld", mintime=5000, combine_avg=True)


    # dir = "/home/victor/Desktop/PhD/src/flatPlateRe3000inc/directLinearSolver/blowingSuction/"
    # combineAvg(dir, chkfile_ending_pattern="_stress.fld", mintime=5000, combine_avg=True)

    ##### multiple cases at once

    # dir_local = "/home/victor/Desktop/PhD/src/incGapRe3000/directLinearSolver/blowingSuction/"
    #
    # # get all directories in dir_local such that dir/data/ does not exist
    #
    # directories = [
    #     # d for d in os.listdir(dir_local) if os.path.isdir(os.path.join(dir_local, d)) and not os.path.exists(os.path.join(dir_local, d, "data"))
    #     # "d0.5_w50",
    #     # "d0.75_w62",
    #     # "d1.5_w70",
    #     # "d2.25_w16",
    #     # "d2.5_w31",
    #     # "d3.5_w21",
    #     # "d1.5_w90",
    #     # "d3.25_w20",
    #     # "d4_w18",
    #     # "d1_w83",
    #     # "d1_w90",
    #     # "d1_w71",
    #     # "d1_w58",
    #     # "d1_w39",
    #     # "d1_w31",
    #     # "d2_w36",
    #     # "d2_w28",
    # ]
    #
    # for case in directories:
    #     dir = os.path.join(
    #         dir_local,
    #         case,
    #     )
    #
    #     combineAvg(
    #         dir,
    #         chkfile_ending_pattern="_stress.fld",
    #         mintime=5000,
    #         combine_avg=True,
    #     )
