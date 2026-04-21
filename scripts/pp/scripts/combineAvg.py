import os
from pp.ssh.fieldconvert import combineAvg


if __name__ == "__main__":
    case = "d4_w17/omega0.08"
    dir = f"/home/victor/Desktop/PhD/src/incGapRe1000/directLinearSolver/omegaBlowSuct/{case}/"
    
    # dir = "/home/victor/Desktop/PhD/src/flatPlateRe1000inc/directLinearSolver/omegaBlowSuct/omega0.65/"


    dir_local = os.path.join("/rds/general/user/vb824/home/Desktop/PhD/runs", dir)


    combineAvg(dir_local, chkfile_ending_pattern="_stress.fld", mintime=5000, combine_avg=True)
