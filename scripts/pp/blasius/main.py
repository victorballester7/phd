import numpy as np
from pp.parameters import Parameters 
import fitting as ft
import plot as pt
import solve_bvp as sbvp


def main():
    param = Parameters(
        eta_int_max=11.3,
        n=401,
        p=11,
        eta_interpolation_max=11,
        x_inflow=-100,
        incNS=False,
        uinf=1.0,
        rhoinf=1.0,
        re_deltaStar=1000,
        mainf=0.4,
        pr=0.72,
        tinf=1.0,
        tref_inf=225.0,  # Tref in Kelvin for the Sutherland law
        twall_dimensionless=1.0,  # = Tinf
        gamma=1.4,
        adiabatic=False,  # False => isothermal wall (Dirichlet BC for T), True => adiabatic wall (Neumann BC for T)
    )

    x = np.linspace(0, param.eta_int_max, param.n)  # Initial mesh

    # Solve the BVP
    y = sbvp.solve_BVP(x, param).sol(x)

    param.setxy_interp(x, y)

    # compute BL parameters
    param.computeBLParams(x, y)
    print(param)

    fields_fits = ft.fitting(x, y, param)

    pt.plot(x, y, fields_fits, param)

    print("FOR SOME REASON THE ERROR ON RHO CHANGES FROM RUN TO RUN (SOMETIMES THE ONE IN ROHU AS WELL). JUST RUN THE PROGRAM UNTIL YOU GET A SMALL ERROR. IT IS A PROBLEM WITH THE CURVE FITTING, NOT WITH THE BVP SOLVER.")


if __name__ == "__main__":  # Physical parameters
    main()
