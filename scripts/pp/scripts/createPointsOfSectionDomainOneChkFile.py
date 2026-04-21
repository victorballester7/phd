import os
import numpy as np
from pp.colors import colors
from pp.ssh.fieldconvert import fld2datapts



def gen_points_file(ymin: float, ymax: float, n: int, dir_local: str, output_file: str):
    pts_path = os.path.join(dir_local, output_file + ".pts")
    
    base_dir = os.path.dirname(pts_path)

    os.makedirs(base_dir, exist_ok=True)

    x_locations_first_minus_eps = np.arange(-70, 100, 1)
    x_locations_first_plus_eps = x_locations_first_minus_eps + 1
    x_locations_first = np.sort(np.unique(np.concatenate((x_locations_first_minus_eps, x_locations_first_plus_eps))))

    x_locations_last_minus_eps = np.arange(150, 1000, 50)
    x_locations_last_plus_eps = x_locations_last_minus_eps + 1
    x_locations_last = np.sort(np.unique(np.concatenate((x_locations_last_minus_eps, x_locations_last_plus_eps))))
    
    if (len(x_locations_first) + len(x_locations_last)) % 2 != 0:
        x_locations_first = x_locations_first[:-1]  # Remove the last point to make it even

    x_locations = np.unique(np.concatenate((x_locations_first, x_locations_last)))
    y_locations = ymin + (ymax - ymin) * (np.linspace(0, 1, n) ** 2)

    X, Y = np.meshgrid(x_locations, y_locations, indexing='ij')
    xy_mesh = np.column_stack((X.ravel(), Y.ravel()))

    print(colors.OKBLUE + f"Generating PTS file with {len(x_locations)} x-locations and {n} y-locations ({len(x_locations)*n} total points)..." + colors.ENDC)

    np.savetxt(pts_path, xy_mesh, fmt="%.6f", header='<?xml version="1.0" encoding="utf-8" ?>\n<NEKTAR>\n<POINTS DIM="2" FIELDS="">', footer='</POINTS>\n</NEKTAR>\n',comments='')

    print(colors.OKGREEN + f"PTS file generated at: {pts_path}" + colors.ENDC)


if __name__ == "__main__":
    case = "d2_w24/omega0.16"
    dir = f"incGapRe1000/directLinearSolver/omegaBlowSuct/{case}/"
    n = 600

    fld_remote = "mesh_avg.fld"
    dir_local = os.path.join("/home/victor/Desktop/PhD/src", dir)

    output_file = f"data/pointsavg_n{n}"


    gen_points_file(ymin=0, ymax=150, n=n, dir_local=dir_local, output_file=output_file)
    fld2datapts(dir_local, fld_remote, output_file)
