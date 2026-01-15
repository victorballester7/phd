from typing import Tuple
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as la
from pp.parameters import Parameters
from scipy.integrate import cumulative_trapezoid
from pp.fileManagement import readFieldsBySection
from pp.filterData import getRMS
import os

def plot(x: np.ndarray, y: np.ndarray, fields_fits: Tuple, param: Parameters) -> None:
    if param.incNS:
        print(
            "\nUsing "
            + "\033[92m"
            + "INCOMPRESSIBLE"
            + "\033[0m"
            + " Navier-Stokes equations\n"
        )

        plot_incNS(x, y, fields_fits, param)
    else:
        print(
            "\nUsing "
            + "\033[92m"
            + "COMPRESSIBLE"
            + "\033[0m"
            + " Navier-Stokes equations\n"
        )

        plot_comNS(x, y, fields_fits, param)


def plot_incNS(
    x: np.ndarray, y: np.ndarray, fields_fits: Tuple, param: Parameters
) -> None:
    # plot the solutions
    ufit, vfit = fields_fits
    re_x = param.getRe_x(param.re_deltaStar)

    # approximate solutions
    u_apx = param.uinf * ufit
    u_apx = np.append(u_apx, param.uinf * np.ones(len(x) - len(param.x_interp)))
    v_apx = 0.5 * param.uinf / np.sqrt(re_x) * vfit
    v_apx = np.append(
        v_apx,
        0.5
        * param.uinf
        / np.sqrt(re_x)
        * param.deltaStar
        * np.ones(len(x) - len(param.x_interp)),
    )

    # derivatives of approximate solutions
    du_apx = np.gradient(u_apx, x, edge_order=2)
    dv_apx = np.gradient(v_apx, x, edge_order=2)

    # exact solutions
    u = param.uinf * y[1]
    v = 0.5 * param.uinf / np.sqrt(re_x) * (x * y[1] - y[0])

    # derivatives of exact solutions
    du = param.uinf * y[2]
    dv = 0.5 * param.uinf / np.sqrt(re_x) * x * y[2]

    # Create the figure and the primary axis
    _, ax1 = plt.subplots()
    ax2 = ax1.twiny()
    plt.grid()
    _, axdiff1 = plt.subplots()
    axdiff2 = axdiff1.twiny()
    plt.grid()

    # Plot u with the first axis
    ax1.plot(u, x, label=r"u exact")
    ax1.axvline(param.uinf, color="black", linestyle="--")
    ax1.set_xlabel(r"$U_\infty f'(x)$")  # same as for v, because we align the axis

    ax2.plot(v, x, label=r"v exact", linestyle="dashdot")

    # Plot diff u in the second figure
    axdiff1.plot(du, x, label="u' exact")
    axdiff1.axvline(0, color="black", linestyle="--")
    axdiff1.set_xlabel("u'")

    axdiff2.plot(dv, x, label="v' exact", linestyle="dashdot")
    axdiff2.set_xlabel("v'")

    uNormLoo, duNormLoo = la.norm(u, np.inf), la.norm(du, np.inf)
    vNormLoo, dvNormLoo = la.norm(v, np.inf), la.norm(dv, np.inf)
    uNormC1 = uNormLoo + duNormLoo
    vNormC1 = vNormLoo + dvNormLoo
    print("||u||_C1 = ||u||_oo + ||u'||_oo")
    print("||u||_oo  = ", uNormLoo)
    print("||u'||_oo = ", duNormLoo)
    print("||u||_C1  = ", uNormC1)
    print()
    print("||v||_oo  = ", vNormLoo)
    print("||v'||_oo = ", dvNormLoo)
    print("||v||_C1  = ", vNormC1)

    erruLoo = la.norm(u - u_apx, np.inf)
    errvLoo = la.norm(v - v_apx, np.inf)

    errduLoo = la.norm(du - du_apx, np.inf)
    errdvLoo = la.norm(dv - dv_apx, np.inf)

    argmax = np.argmax(np.abs(du - du_apx))
    print("max of error du: ", x[argmax])
    argmax = np.argmax(np.abs(dv - dv_apx))
    print("max of error dv: ", x[argmax])

    print("---- u ----")
    print("||u - u_aprox||_oo / ||u||_oo    = ", erruLoo / uNormLoo)
    print("||u' - u_aprox'||_oo / ||u'||_oo = ", errduLoo / duNormLoo)
    print(
        "||u - u_aprox||_C1 / ||u||_C1 (error below is always greater than this)  = ",
        (erruLoo + errduLoo) / uNormC1,
    )
    print(
        "max(||u - u_aprox||_oo / ||u||_oo, ||u' - u_aprox'||_oo / ||u'||_oo)     = ",
        np.fmax(erruLoo / uNormLoo, errduLoo / duNormLoo),
    )
    print()
    print("---- v ----")
    print("||v - v_aprox||_oo / ||v||_oo    = ", errvLoo / vNormLoo)
    print("||v' - v_aprox'||_oo / ||v'||_oo = ", errdvLoo / dvNormLoo)
    print(
        "||v - v_aprox||_C1 / ||v||_C1 (error below is always greater than this)    = ",
        (errvLoo + errdvLoo) / vNormC1,
    )
    print(
        "max(||v - v_aprox||_oo / ||v||_oo, ||v' - v_aprox'||_oo / ||v'||_oo)     = ",
        np.fmax(errvLoo / vNormLoo, errdvLoo / dvNormLoo),
    )

    ax1.plot(u_apx, x, label="u fit")
    ax2.plot(v_apx, x, label="v fit", linestyle="dashdot")

    axdiff1.plot(du_apx, x, label="u' fit")
    axdiff2.plot(dv_apx, x, label="v' fit", linestyle="dashdot")

    eps = 0.05
    x1_min = 0
    x1_max = u[-1]
    x1_min -= eps * (x1_max - x1_min)
    x1_max += eps * (x1_max - x1_min)
    ax1.set_xlim(x1_min, x1_max)

    x2_min = 0
    x2_max = v[-1]
    x2_min -= eps * (x2_max - x2_min)
    x2_max += eps * (x2_max - x2_min)
    ax2.set_xlim(x2_min, x2_max)

    ax1.legend(loc="upper left")
    ax2.legend(loc="upper right")
    axdiff1.legend(loc="upper left")
    axdiff2.legend(loc="upper right")
    plt.show()


def plot_comNS(
    x: np.ndarray, y: np.ndarray, fields_fits: Tuple, param: Parameters
) -> None:
    rho_fit, rhou_fit, rhov_fit, e1_fit, e2_fit = fields_fits

    # plot the solutions
    re_x = param.getRe_x(param.re_deltaStar)
    rho_apx = param.rhoinf * rho_fit
    rho_apx = np.append(rho_apx, param.rhoinf * np.ones(len(x) - len(param.x_interp)))
    rhou_apx = param.rhoinf * param.uinf * rhou_fit
    rhou_apx = np.append(
        rhou_apx, param.rhoinf * param.uinf * np.ones(len(x) - len(param.x_interp))
    )
    rhov_apx = 0.5 * param.rhoinf * param.uinf / np.sqrt(re_x) * rhov_fit
    rhov_apx = np.append(
        rhov_apx,
        0.5
        * param.rhoinf
        * param.uinf
        / np.sqrt(re_x)
        * param.deltaStar
        * np.ones(len(x) - len(param.x_interp)),
    )
    e_apx = (
        param.rhoinf
        * param.uinf**2
        * (
            0.5 * e1_fit
            + 0.125 / re_x * e2_fit
            + 1.0 / ((param.gamma - 1) * param.gamma * param.mainf**2)
        )
    )
    e_apx = np.append(
        e_apx,
        param.rhoinf
        * param.uinf**2
        * (
            0.5
            + 0.125 / re_x * param.deltaStar**2
            + 1.0 / ((param.gamma - 1) * param.gamma * param.mainf**2)
        )
        * np.ones(len(x) - len(param.x_interp)),
    )

    # exact solutions
    rho = param.rhoinf / y[3]
    rhou = param.rhoinf * param.uinf * y[1] / y[3]
    rhov = (
        0.5
        * param.rhoinf
        * param.uinf
        / np.sqrt(re_x)
        * (y[1] / y[3] * cumulative_trapezoid(y[3], x, initial=0) - y[0])
    )
    e = (
        param.rhoinf
        * param.uinf**2
        * (
            0.5 * y[1] ** 2 / y[3]
            + 0.125
            / re_x
            * (
                y[1] / np.sqrt(y[3]) * cumulative_trapezoid(y[3], x, initial=0)
                - np.sqrt(y[3]) * y[0]
            )
            ** 2
            + 1.0 / ((param.gamma - 1) * param.gamma * param.mainf**2)
        )
    )

    # === Norms and errors ===
    real = [rho, rhou, rhov, e]
    approx = [rho_apx, rhou_apx, rhov_apx, e_apx]
    labels = ["ρ", "ρu", "ρv", "E"]

    deriv_real = []
    deriv_approx = []

    for r, a, lab in zip(real, approx, labels):
        # norms L∞
        normLoo = la.norm(r, np.inf)
        d_r = np.gradient(r, x, edge_order=2)
        normdLoo = la.norm(d_r, np.inf)
        normC1 = normLoo + normdLoo

        # errors
        errLoo = la.norm(r - a, np.inf)
        d_a = np.gradient(a, x, edge_order=2)
        errdLoo = la.norm(d_r - d_a, np.inf)

        argmax = np.argmax(np.abs(d_r - d_a))
        print(f"\n--- {lab} ---")
        print(f"||{lab}||_oo   = {normLoo}")
        print(f"||{lab}'||_oo  = {normdLoo}")
        print(f"||{lab}||_C1   = {normC1}")
        print(f"||{lab} - {lab}_apx||_oo / ||{lab}||_oo    = {errLoo / normLoo}")
        print(f"||{lab}' - {lab}_approx'||_oo / ||{lab}'||_oo = {errdLoo / normdLoo}")
        print(
            f"||{lab} - {lab}_approx||_C1 / ||{lab}||_C1 (error below is always greater than this)  = {(errLoo + errdLoo) / normC1}"
        )
        print(
            f"max(||{lab} - {lab}_approx||_oo / ||{lab}||_oo, ||{lab}' - {lab}_approx'||_oo / ||{lab}'||_oo)     = {np.fmax(errLoo / normLoo, errdLoo / normdLoo)}"
        )
        print(f"max of error {lab}' at x = {x[argmax]}")

        deriv_real.append(d_r)
        deriv_approx.append(d_a)

    # normalize energy for plotting
    e = e / e[-1]
    e_apx = e_apx / e_apx[-1]

    # === Plots: fields ===
    _, ax1 = plt.subplots()
    ax1.plot(rho, x, label="ρ exact", color="orange")
    ax1.plot(rho_apx, x, label="ρ fit", color="orange", linestyle="dotted")
    ax1.plot(rhou, x, label="ρu exact", color="blue")
    ax1.plot(rhou_apx, x, label="ρu fit", color="blue", linestyle="dotted")
    ax1.plot(e, x, label="E / lim E exact", color="purple")
    ax1.plot(e_apx, x, label="E / lim E fit", color="purple", linestyle="dotted")
    ax1.axvline(param.uinf, color="black", linestyle="--")

    ax2 = ax1.twiny()
    ax2.plot(rhov, x, label="ρv exact", color="green")
    ax2.plot(rhov_apx, x, label="ρv fit", color="green", linestyle="dotted")

    ##################### tmp code to plot the section of the doamin from DNS data ###############

        # basePath = "../../../src/flatSurfaceRe1000Ma0.6ComNS/dns/"
    # basePath = os.path.join(os.path.dirname(os.path.abspath(__file__)), basePath)
    # chkfile = "13"
    # n = 800
    # dataFile = os.path.join(basePath, "data", f"points{chkfile}_all_n{n}.dat")

    # xvals, yvals, data = readFieldsBySection(dataFile)

    # rms = getRMS(data)


    # this is temporary, comment to plot reynold stresses
   #  rhoo = data[:, :, 2]
   #  u = data[:, :, 3]
   #  v = data[:, :, 4]
   #  rms = u
   # 
   #  etaMax = 10.98925
   #  yMax = 5.202553147076392 

   #  # def eta(y):
   #  #     if y < yMax:
   #  #         return 2.120192718897813 * (y + a1_y2eta * np.tanh(b1_y2eta * y) + a2_y2eta * np.tanh(b2_y2eta * y) ** 2)
   #  #     else:
   #  #         return etaMax

   #  a0_y2eta = 2.120192718897813
   #  a1_y2eta = 0.007174212003742566
   #  b1_y2eta = 0.35787752278060175
   #  a2_y2eta = -0.0480154451948319
   #  b2_y2eta = 0.9143119449738151
   #  def eta(y):
   #      if y < yMax:
   #          return a0_y2eta * y + a1_y2eta *np.tanh(b1_y2eta * y) * a2_y2eta * np.tanh(b2_y2eta * y)**2
   #      else:
   #          return etaMax

    # def rhofit(y, a_coeffs, b_coeffs):
    #     e = eta(y)
    #     num = sum(a * e**(i+1) for i, a in enumerate(a_coeffs))
    #     den = 1 + sum(b * e**(i+1) for i, b in enumerate(b_coeffs))
    #     return 1 + num / den

    # def ufit(y, a_coeffs, b_coeffs):
    #     e = eta(y)
    #     num = sum(a * e**(i+1) for i, a in enumerate(a_coeffs))
    #     den = 1 + sum(b * e**(i+1) for i, b in enumerate(b_coeffs))
    #     return num / den

    # rho = np.array([rhofit(y, a_rho, b_rho) for y in yvals])
    # rhou = np.array([ufit(y, a_rhou, b_rhou) for y in yvals])
    # 
    # ax.plot(rho, yvals, "k--", label="u profile")
    # ax.plot(eta(yvals), yvals, "r--", label="eta profile")
    
    # yvals = yvals / yMax * etaMax

    # indx_eta = np.argmin(np.abs(yvals - 12))
    # indx_eta = -1
    # yvals = yvals[:indx_eta]
    # rms = rms[:, :indx_eta]
    # for i,y in enumerate(yvals):
    #     yvals[i] = eta(y)


    # # Plot the points
    # for i in range(len(xvals)):
    #     ax1.plot(rms[i], yvals, "-", markersize=2, label=f"x = {xvals[i]}")


    ##############################################################################################

    eps = 0.05
    x1_min = 0
    x1_max = rhou[-1]
    x1_min -= eps * (x1_max - x1_min)
    x1_max += eps * (x1_max - x1_min)
    ax1.set_xlim(x1_min, x1_max)

    x2_min = 0
    x2_max = rhov[-1]
    x2_min -= eps * (x2_max - x2_min)
    x2_max += eps * (x2_max - x2_min)
    ax2.set_xlim(x2_min, x2_max)

    ax1.legend(loc="upper left")
    ax2.legend(loc="upper right")
    plt.grid()

    # === Plots: derivatives ===
    _, axdiff1 = plt.subplots()
    axdiff2 = axdiff1.twiny()

    for d_r, d_a, lab, color in zip(
        deriv_real,
        deriv_approx,
        labels,
        ["orange", "blue", "green", "purple"],
    ):
        if lab != "ρv":
            axdiff1.plot(d_r, x, label=f"{lab}' exact", color=color)
            axdiff1.plot(d_a, x, label=f"{lab}' fit", color=color, linestyle="dotted")
        else:
            axdiff2.plot(d_r, x, label=f"{lab}' exact", color=color)
            axdiff2.plot(d_a, x, label=f"{lab}' fit", color=color, linestyle="dotted")

    axdiff1.axvline(0, color="black", linestyle="--")


    eps = 0.15
    x1_min = 0
    x1_max = np.max(deriv_approx[1])
    x1_min -= eps * (x1_max - x1_min)
    x1_max += eps * (x1_max - x1_min)
    axdiff1.set_xlim(x1_min, x1_max)

    x2_min = 0
    x2_max = np.max(deriv_approx[2])
    x2_min -= eps * (x2_max - x2_min)
    x2_max += eps * (x2_max - x2_min)
    axdiff2.set_xlim(x2_min, x2_max)
    axdiff1.legend(loc="upper left")
    axdiff2.legend(loc="upper right")
    plt.grid()
    plt.show()
