import os
import numpy as np
from pp.colors import colors
from pp.ssh.fieldconvert import fld2datapts_batch
from pp.fileManagement import extract_depth_width


def gen_points_file(
    x: float, ymin: float, ymax: float, n: int, dir_local: str, output_file: str
):
    pts_path = os.path.join(dir_local, output_file + ".pts")

    base_dir = os.path.dirname(pts_path)

    os.makedirs(base_dir, exist_ok=True)

    # x_locations = np.unique(np.concatenate((x_locations_first, x_locations_last)))
    y_locations = ymin + (ymax - ymin) * (np.linspace(0, 1, n) ** 2)

    X, Y = np.meshgrid(x, y_locations, indexing="ij")
    xy_mesh = np.column_stack((X.ravel(), Y.ravel()))

    print(
        colors.OKBLUE
        + f"Generating PTS file with 1 x-locations and {n} y-locations)..."
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
    print(f"Width: {w}, Depth: {d}")
    n = 600

    # x_loc = np.arange(-70, 1000, 10) + 3.5 # 3.5 just not to coincide with the edges of the gap geoemtry
    if d != 0:
        x_loc1 = np.arange(-70, -10, 20) + 3.5
        x_loc2 = (
            np.arange(-10, w + 40, 1) + 1.5
        )  # 3.5 just not to coincide with the edges of the gap geoemtry
        x_loc3 = np.arange(w + 40, 1000, 20) + 3.5
        x_loc = np.unique(np.concatenate((x_loc1, x_loc2, x_loc3)))
        # x_loc = x_loc2
    else:
        x_loc = (
            np.arange(-70, 1000, 20) - 3.5
        )  # 3.5 just not to coincide with the edges of the gap geoemtry
    fld_remote = "baseflow.fld"
    ymax = 150

    output_files = []

    for x in x_loc:
        output_file = f"data/xSections/points_n{n}_x{x}"
        if (x > 0 and x < w) or (x > 0 and w == 0):
            ymin = -d
        else:
            ymin = 0

        gen_points_file(
            x, ymin, ymax, n=n, dir_local=dir_local, output_file=output_file
        )
        output_files.append(output_file)

    fld2datapts_batch(dir_local, fld_remote, output_files)


if __name__ == "__main__":
    # dir_local = "/home/victor/Desktop/PhD/src/flatPlateRe1000inc/directLinearSolver/blowingSuction/wgnInsideDomainDivFree/"
    # dir_local = "/home/victor/Desktop/PhD/src/incGapRe1000/directLinearSolver/blowingSuctionCoarserMesh/"

    # # case = ""
    # case = "d2.25_w31"
    # create_points_file(dir_local, case=case)

    ##### multiple cases at once

    dir_local = (
        # "/home/victor/Desktop/PhD/src/bfsRe1000inc/directLinearSolver/blowingSuction/"
        # "/home/victor/Desktop/PhD/src/incGapRe1000/directLinearSolver/blowingSuction/"
        "/home/victor/Desktop/PhD/src/incGapRe800/directLinearSolver/blowingSuction/"
    )

    # get all directories in dir_local such that dir/data/ does not exist

    directories = [
        # "d0.5_w90",
        "d4_w12",
        # for d in os.listdir(dir_local)
        # if os.path.isdir(os.path.join(dir_local, d))
        # and not os.path.exists(os.path.join(dir_local, d, "data"))
        ]

    for case in directories:
        create_points_file(dir_local, case=case)
