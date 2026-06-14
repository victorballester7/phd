import os
import numpy as np
from pp.colors import colors
from pp.ssh.fieldconvert import fld2datapts


def gen_points_file(ymin: float, ymax: float, n: int, dir_local: str, output_file: str):
    pts_path = os.path.join(dir_local, output_file + ".pts")

    base_dir = os.path.dirname(pts_path)

    os.makedirs(base_dir, exist_ok=True)

    x_locations = np.arange(-70, 1000, 10) + 3.5
   
    # x_locations = np.unique(np.concatenate((x_locations_first, x_locations_last)))
    y_locations = ymin + (ymax - ymin) * (np.linspace(0, 1, n) ** 2)

    X, Y = np.meshgrid(x_locations, y_locations, indexing="ij")
    xy_mesh = np.column_stack((X.ravel(), Y.ravel()))

    print(
        colors.OKBLUE
        + f"Generating PTS file with {len(x_locations)} x-locations and {n} y-locations ({len(x_locations) * n} total points)..."
        + colors.ENDC
    )

    np.savetxt(
        pts_path,
        xy_mesh,
        fmt="%.6f",
        header='<?xml version="1.0" encoding="utf-8" ?>\n<NEKTAR>\n<POINTS DIM="2" FIELDS="">',
        footer="</POINTS>\n</NEKTAR>\n",
        comments="",
    )

    print(colors.OKGREEN + f"PTS file generated at: {pts_path}" + colors.ENDC)


def create_points_file(base_dir: str, case: str):
    dir_local = base_dir + case
    n = 600
    output_file = f"data/pointsavg_n{n}"
    fld_remote = "mesh_avg.fld"
    gen_points_file(
        ymin=-4, ymax=150, n=n, dir_local=dir_local, output_file=output_file
    )
    fld2datapts(dir_local, fld_remote, output_file)


if __name__ == "__main__":
    # case = "d2_w24/omega0.08"
    # dir_local = f"/home/victor/Desktop/PhD/src/incGapRe1000/directLinearSolver/omegaBlowSuct/{case}"
    # dir_local = "/home/victor/Desktop/PhD/src/flatPlateRe1000inc/directLinearSolver/omegaBlowSuct/omega0.04"


    # dir = "/home/victor/Desktop/PhD/src/bfsRe1000inc/directLinearSolver/blowingSuction/"
    # case = "d1.5"
    # create_points_file(dir, case=case)



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

    # directories = [
    #     # "d1.5_w200", "d1_w200", "d0.5_w200", 
    #     d for d in os.listdir(dir_local) if os.path.isdir(os.path.join(dir_local, d)) and not os.path.exists(os.path.join(dir_local, d, "data"))
    # ]

    for case in directories:
        create_points_file(dir_local, case=case)
