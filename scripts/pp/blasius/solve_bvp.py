import numpy as np
from scipy.integrate import _bvp, solve_bvp
from scipy.interpolate import PPoly

from pp.parameters import Parameters


def solve_BVP(x: np.ndarray, param: Parameters):
    if param.incNS:
        return solve_BVP_incNS(x)
    else:
        return solve_BVP_comNS(x, param)


def solve_BVP_incNS(x: np.ndarray):
    # Define the system of ODEs
    def odes(x, y):
        f0, f1, f2 = y
        df0_dx = f1
        df1_dx = f2
        df2_dx = -0.5 * f0 * f2
        return np.vstack((df0_dx, df1_dx, df2_dx))

    # Define the boundary conditions
    def bc(ya, yb):
        return np.array(
            [
                ya[0],  # y1(0) = f(0) = 0
                ya[1],  # y2(0) = f'(0) = 0
                yb[1] - 1,  # y2(inf) = f'(inf) = 1
            ]
        )

    y_guess = np.zeros((3, x.size))
    eta_max = x[-1]
    y_guess[1, :] = x / eta_max  # Approximate guess for f'(x)

    # Solve the BVP
    solution = solve_bvp(odes, bc, x, y_guess, tol=1e-10, max_nodes=10000)


    if solution.success:
        print("BVP solved successfully!")
    else:
        print("BVP solution failed.")

    return solution


def solve_BVP_comNS(x: np.ndarray, param: Parameters):
    # ref values for Sutherland's law:
    # mu = muRef * (T/TRef)**1.5 * (TRef + S)/(T + S)

    # dimensionless mu (with T dimensionless)
    # mu = T**1.5 * (1 + S/Tref) / (T + S/Tref) 
    # where T is dimensionless and Tref is in Kelvin (I think Tref as dimensional Tinf)

    # g := 1/rho (rho dimensionless, e.g. rho := tilde{rho} := rho/rhoInf)
    # since p_x = 0 (assumption for self-similarity solution), we have that in the BL p=const ( because from the v-component of the momentum eq we have p_y = 0), thus from perfect gas law (p = rho R T) we have that rhoInf/rho = T/Tinf ==> g = T (dimensionless)
    # C := mu rho = mu / g = g**0.5 *(1 + S/TInf) / (g + S/TInf)

    # muRef = 1.716e-5
    # TRef = 273.15
    S = 110.0 # same as in Nektar++

    # Define the system of ODEs
    def odes(x, y):
        f0, f1, f2, g0, g1 = y

        C = g0**0.5 * (1 + S / param.tref_inf) / (g0 + S / param.tref_inf)
        dC = (0.5 / g0 - 1.0 / (g0 + S / param.tref_inf)) * g1 * C

        df0_dx = f1
        df1_dx = f2
        df2_dx = -f2 * (dC + f0 / 2) / C
        dg0_dx = g1
        dg1_dx = (
            -g1 * (dC + param.pr * f0 / 2) / C
            - (param.gamma - 1) * param.pr * param.mainf**2 * f2**2
        )
        return np.vstack((df0_dx, df1_dx, df2_dx, dg0_dx, dg1_dx))

    # Define the boundary conditions
    def bc(ya, yb):
        bcs = np.array(
            [
                ya[0],  # y1(0) = f(0) = 0
                ya[1],  # y2(0) = f'(0) = 0
                yb[1] - 1,  # y2(inf) = f'(inf) = 1
                yb[3] - 1,  # g(inf) = 1 (rho(inf) = rhoInf)
                # and exactly one of the two below (first is for isothermal - dirichlet BC, second for adiabatic - Neumann BC)
                # ya[3] - 1,  # g(0) = 1 (if T_wall = TInf)
                # ya[4],  # g'(0) = 0 (no heat flux at the wall)
                # it is decided based on param.adiabatic
            ]
        )

        if param.adiabatic:
            bcs = np.append(bcs, ya[4])
        else:  # isothermal
            bcs = np.append(
                bcs, ya[3] - param.twall_dimensionless
            )  # g(0) = T_wall / TInf
        return bcs

    y_guess = np.ones((5, x.size)) / 2
    y_guess[0] = x**2 / (2 * x[-1])  # f(eta)
    y_guess[1] = x / x[-1]  # f'
    y_guess[2] = 1 / x[-1]  # f''
    y_guess[3] = 1 - (1 - param.twall_dimensionless) * np.exp(-x)  # T(eta)
    y_guess[4] = (1 - param.twall_dimensionless) * np.exp(-x)  # T'(eta)

    # Solve the BVP
    solution = solve_bvp(odes, bc, x, y_guess, tol=1e-10, max_nodes=10000)

    if solution.success:
        print("BVP solved successfully!")
    else:
        print("BVP solution failed.")

    return solution
