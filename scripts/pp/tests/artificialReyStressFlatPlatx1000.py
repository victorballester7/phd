import numpy as np
from pp.fileManagement import readFieldsBySection
from pp.colors import colors



def main():
    dir = "/home/victor/Desktop/PhD/src/flatPlateRe1000inc/directLinearSolver/blowingSuction/data"
    dataFile_old = dir + "/pointsavg_n600.dat.bak"
    dataFile_new = dir + "/pointsavg_n600.dat"

    try:
        x, y, data = readFieldsBySection(dataFile_old)
    except Exception as e:
        print(colors.FAIL + f"Error reading data from {dataFile_old}: {e}" + colors.ENDC)
        return np.array([]), np.array([])

    nx_850 = 2.64061815
    nx_900 = 2.72308841
    nx_950 = 2.80233110

    # n_1000_approx = 2*nx_950 - nx_900
    n_1000_approx = 5/2. * nx_950 - 2 * nx_900 + 1/2. * nx_850
    alpha = np.exp(n_1000_approx - nx_950)
    beta = alpha
    
    # line by line copy old file to new file, except for when x = 1000, in which case multiply the data by beta
    with open(dataFile_old, "r") as f_old, open(dataFile_new, "w") as f_new:
        for i, line in enumerate(f_old):
            if i < 3:
                f_new.write(line)
                continue
            parts = line.split()
            x_val = float(parts[0])
            if x_val == 1000.0:
                continue
            f_new.write(" ".join(parts) + "\n")

    data_950 = data[x == 950]

    # multiply u, v by sqrt(beta) and uu, vv by beta, and the others by nan to stress that they are not valid
    data_1000 = np.copy(data_950).reshape(data.shape[1],data.shape[2])
        
    print(colors.OKBLUE + f"Scaling data at x=950 to x=1000 using beta={beta:.6f}..." + colors.ENDC)

    data_1000[:, 0] = 1000.0  # x
    data_1000[:, 2] *= np.sqrt(beta)  # u
    data_1000[:, 3] *= np.sqrt(beta)  # v
    data_1000[:, 4] *= np.nan # p
    data_1000[:, 5] *= beta  # uu
    data_1000[:, 6] *= np.nan  # uv
    data_1000[:, 7] *= beta  # vv

    # append data_1000 to the new file
    with open(dataFile_new, "a") as f_new:
        for row in data_1000:
            f_new.write(" ".join(map(str, row)) + "\n")
    
if __name__ == "__main__":
    main()
