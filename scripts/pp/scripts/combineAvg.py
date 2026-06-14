import os
from pp.ssh.fieldconvert import combineAvg


if __name__ == "__main__":
    # case = "d4_w17/omega0.08"
    # dir = f"/home/victor/Desktop/PhD/src/incGapRe1000/directLinearSolver/omegaBlowSuct/{case}/"

    # dir = "/home/victor/Desktop/PhD/src/flatPlateRe1000inc/directLinearSolver/omegaBlowSuct/omega0.02/"

    # dir = "/home/victor/Desktop/PhD/src/bfsRe1000inc/directLinearSolver/blowingSuction/d1.5/"
    # combineAvg(dir, chkfile_ending_pattern="_stress.fld", mintime=5000, combine_avg=True)

    ##### multiple cases at once

    dir_local = "/home/victor/Desktop/PhD/src/incGapRe1000/directLinearSolver/blowingSuctionCoarserMesh/"

    # get all directories in dir_local such that dir/data/ does not exist

    directories = [
        # d for d in os.listdir(dir_local) if os.path.isdir(os.path.join(dir_local, d)) and not os.path.exists(os.path.join(dir_local, d, "data"))
        "d1_w83",
        "d1_w90",
        "d1_w71",
        "d1_w58",
        "d1_w39",
        "d1_w31",
        "d2_w36",
        "d2_w28",
    ]

    for case in directories:
        dir = os.path.join(
            dir_local,
            case,
        )

        combineAvg(
            dir,
            chkfile_ending_pattern="_stress.fld",
            mintime=5000,
            combine_avg=True,
        )
