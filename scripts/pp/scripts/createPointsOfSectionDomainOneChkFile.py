import os
import numpy as np
from pp.colors import colors
from pp.ssh.fieldconvert import fld2datapts
from pp.fileManagement import extract_depth_width


def gen_points_file(
    ymin: float, ymax: float, n: int, width: float, dir_local: str, output_file: str
):
    pts_path = os.path.join(dir_local, output_file + ".pts")

    base_dir = os.path.dirname(pts_path)

    os.makedirs(base_dir, exist_ok=True)

    x_locations = np.arange(width, width + 950, 2) + 1.5
    # x_loc1 = np.arange(-70, -10, 10) + 3.5
    # x_loc2 = (
    #     np.arange(-10, width + 100, 2) + 1.5
    # )  # 3.5 just not to coincide with the edges of the gap geoemtry
    # x_loc3 = np.arange(width + 100, 1000, 10) + 3.5
    # x_locations = np.unique(np.concatenate((x_loc1, x_loc2, x_loc3)))
    # x_loc = x_loc2

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
    d, w = extract_depth_width(case)
    n = 600
    # output_file = f"data/pointsPO_n{n}"
    # fld_remote = "mesh_po.fld"
    output_file = f"data/pointsavg_n{n}"
    fld_remote = "mesh_avg.fld"
    gen_points_file(
        # ymin=0, ymax=150, n=n, width=w, dir_local=dir_local, output_file=output_file
        ymin=-d, ymax=150, n=n, width=w, dir_local=dir_local, output_file=output_file
    )
    fld2datapts(dir_local, fld_remote, output_file)


if __name__ == "__main__":
    # case = "d2_w24/omega0.08"
    # dir = "/home/victor/Desktop/PhD/src/incGapRe3000/directLinearSolver/blowingSuction/"
    # dir_local = "/home/victor/Desktop/PhD/src/flatPlateRe1000inc/directLinearSolver/omegaBlowSuct/omega0.04"

    # dir = "/home/victor/Desktop/PhD/src/bfsRe1000inc/directLinearSolver/blowingSuction/"
    # create_points_file(dir, case="d0.25_w40")

    # dir = "/home/victor/Desktop/PhD/src/flatPlateRe3000inc/directLinearSolver/blowingSuction/"
    # create_points_file(dir, case="")

    ##### multiple cases at once

    # dir_local = "/home/victor/Desktop/PhD/src/quasi3DGapRe1000inc/baseflow/dns/"
    # dir_local = "/home/victor/Desktop/PhD/src/incGapRe1000/baseflow/dns/"
    dir_local = "/home/victor/Desktop/PhD/src/incGapRe1000/directLinearSolver/blowingSuction/"
    #
    # # get all directories in dir_local such that dir/data/ does not exist
    #
    directories = [
        # d
        # for d in os.listdir(dir_local)
        # if os.path.isdir(os.path.join(dir_local, d))
        # and not os.path.exists(os.path.join(dir_local, d, "data"))
        # "d1.75_w40_new",
        # "d1.5_w60_Lz5_k32",
        "d1.5_w64",
        # "d1.5_w63_Lz6.0_k32",
        # "d1.5_w65_Lz6.0_k32",
        # "d1.5_w70_Lz6.0_k32",
        # "d1.5_w80_Lz6.0_k32"
    ]

    for case in directories:
        # print(f"Creating points file for case: {case}")
        print(
            colors.OKBLUE
            + f"Creating points file for case {np.where(np.array(directories) == case)[0][0] + 1}/{len(directories)}: {case}"
            + colors.ENDC
        )
        create_points_file(dir_local, case=case)
