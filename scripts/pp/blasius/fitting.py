from typing import Tuple
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import cumulative_trapezoid
import scipy.linalg as la
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from pp.parameters import Parameters 


def fitting(x: np.ndarray, y: np.ndarray, param: Parameters) -> Tuple[np.ndarray, ...]:
    limit_u = 1.0
    limit_v = param.deltaStar if param.incNS else np.trapezoid(y[3], x) - y[0][-1]

    if param.incNS:
        return fitting_incNS(param.x_interp, param.y_interp, param, [limit_u, limit_v])
    else:
        limit_rho = 1.0
        return fitting_comNS(
            param.x_interp,
            param.y_interp,
            param,
            [
                limit_rho,  # limit_rho
                limit_rho * limit_u,  # limit_rhou
                limit_v,  # limit_rhov
                limit_rho * limit_u**2,  # limit_E1
                limit_v**2,  # limit_E2
            ],
        )


def fitting_incNS(
    x: np.ndarray, y: np.ndarray, param: Parameters, limits: list
) -> Tuple[np.ndarray, np.ndarray]:
    limit_u, limit_v = limits
    pu = np.ones(2 * param.p - 1)  # initial guess for the coefficients
    pv = np.ones(2 * param.p)  # initial guess for the coefficients

    func_tofit_u = y[1]
    func_tofit_v = x * y[1] - y[0]

    ufit_coeffs, _ = curve_fit(
        lambda eta, *params: u_fit_func(eta, param.p, limit_u, *params),
        x,
        func_tofit_u,
        p0=0 * pu,
    )
    ufit = u_fit_func(x, param.p, limit_u, *ufit_coeffs)

    vfit_coeffs, _ = curve_fit(
        lambda eta, *params: v_fit_func(eta, param.p, limit_v, *params),
        x,
        func_tofit_v,
        p0=pv,
    )
    vfit = v_fit_func(x, param.p, limit_v, *vfit_coeffs)

    printCoeffs_incNS([ufit_coeffs, vfit_coeffs], param, [limit_u, limit_v])

    return ufit, vfit


def fitting_comNS(
    x: np.ndarray, y: np.ndarray, param: Parameters, limits: list
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    pu = np.ones(2 * param.p - 1)  # initial guess for the coefficients
    pv = np.ones(2 * param.p)  # initial guess for the coefficients

    limit_rho, limit_rhou, limit_rhov, limit_e1, limit_e2 = limits

    f_tofit_rho = 1 / y[3]
    f_tofit_rhou = y[1] / y[3]
    f_tofit_rhov = y[1] / y[3] * cumulative_trapezoid(y[3], x, initial=0) - y[0]
    # we split the energy equation into two parts:
    f_tofit_e1 = y[1] ** 2 / y[3]
    f_tofit_e2 = (
        y[1] / np.sqrt(y[3]) * cumulative_trapezoid(y[3], x, initial=0)
        - y[0] * np.sqrt(y[3])
    ) ** 2

    rho_fit_coeffs, _ = curve_fit(
        lambda eta, *params: rho_fit_func(eta, param.p, limit_rho, *params),
        x,
        f_tofit_rho,
        p0=pu,
    )
    rhou_fit_coeffs, _ = curve_fit(
        lambda eta, *params: u_fit_func(eta, param.p, limit_rhou, *params),
        x,
        f_tofit_rhou,
        p0=pu,
    )
    rhov_fit_coeffs, _ = curve_fit(
        lambda eta, *params: v_fit_func(eta, param.p, limit_rhov, *params),
        x,
        f_tofit_rhov,
        p0=pv,
    )
    e1_fit_coeffs, _ = curve_fit(
        lambda eta, *params: u_fit_func(eta, param.p, limit_e1, *params),
        x,
        f_tofit_e1,
        p0=pu,
    )
    e2_fit_coeffs, _ = curve_fit(
        lambda eta, *params: v_fit_func(eta, param.p, limit_e2, *params),
        x,
        f_tofit_e2,
        p0=pv,
    )
    y2eta_coeffs = approx_y2eta_comNS(x, y, param)

    printCoeffs_comNS(
        [
            rho_fit_coeffs,
            rhou_fit_coeffs,
            rhov_fit_coeffs,
            e1_fit_coeffs,
            e2_fit_coeffs,
            y2eta_coeffs,
        ],
        param,
        limits,
    )

    rho_fit = rho_fit_func(x, param.p, limit_rho, *rho_fit_coeffs)
    rhou_fit = u_fit_func(x, param.p, limit_rhou, *rhou_fit_coeffs)
    rhov_fit = v_fit_func(x, param.p, limit_rhov, *rhov_fit_coeffs)
    e1_fit = u_fit_func(x, param.p, limit_e1, *e1_fit_coeffs)
    e2_fit = v_fit_func(x, param.p, limit_e2, *e2_fit_coeffs)

    return rho_fit, rhou_fit, rhov_fit, e1_fit, e2_fit


# rational function for f' (both incNS and comNS)
def u_fit_func(
    eta: np.ndarray, p: int, limit_u: float, *params: np.ndarray
) -> np.ndarray:
    num = sum(
        params[i] * eta ** (i + 1) for i in range(p)
    )  # we impose that the function is zero at the origin
    den = (
        1
        + sum(params[p + j] * eta ** (j + 1) for j in range(p - 1))
        + params[p - 1] / limit_u * eta**p
    )  # we impose that the function is 1 at infinity
    return num / den


def v_fit_func(
    eta: np.ndarray, p: int, limit_v: float, *params: np.ndarray
) -> np.ndarray:
    num = sum(
        params[i] * eta ** (i + 2) for i in range(p)
    )  # we impose that the function and its first derivative are zero at the origin
    den = (
        1
        + sum(params[p + j] * eta ** (j + 1) for j in range(p))
        + params[p - 1] / limit_v * eta ** (p + 1)
    )  # we impose that the function is 1 at infinity
    return num / den


def rho_fit_func(
    eta: np.ndarray, p: int, limit_rho: float, *params: np.ndarray
) -> np.ndarray:
    num = 1 + sum(
        params[i] * eta ** (i + 1) for i in range(p)
    )  # we impose that the function is zero at the origin
    den = (
        1
        + sum(params[p + j] * eta ** (j + 1) for j in range(p - 1))
        + params[p - 1] / limit_rho * eta**p
    )  # we impose that the function is 1 at infinity
    result = num / den
    return np.array(result)


def approx_y2eta_comNS(x: np.ndarray, y: np.ndarray, param: Parameters) -> np.ndarray:
    lenParams = 4
    p0 = np.ones(lenParams)  # initial guess for the coefficients

    deltaStar_inflow = np.sqrt(
        1 + param.x_inflow * param.deltaStar**2 / param.re_deltaStar
    )

    # x / sqrt(re_x) = deltaStar_x / 1.72... = deltaStar_x / param.deltaStar
    y_physical_inflow = (
        deltaStar_inflow / param.deltaStar * cumulative_trapezoid(y[3], x, initial=0)
    )

    def int_rho_fit(y, *params):
        alpha, beta, a, b = params
        return y + alpha * np.tanh(beta * y) + a * np.tanh(b * y) ** 2

    f_tofit = cumulative_trapezoid(1 / y[3], y_physical_inflow, initial=0)

    fit_coeffs, _ = curve_fit(int_rho_fit, y_physical_inflow, f_tofit, p0=p0)

    print(f"Approximation of η(x,y) at x_inflow (x = {param.x_inflow}):")
    print(
        f"η(x_inflow,y) = {param.deltaStar / deltaStar_inflow} * (  y + {fit_coeffs[0]} * tanh({fit_coeffs[1]} * y) + {fit_coeffs[2]} * tanh({fit_coeffs[3]} * y)^2 )"
    )

    # plot int_0^y rho dy
    xplot = y_physical_inflow
    yplot_real = cumulative_trapezoid(1 / y[3], y_physical_inflow, initial=0)
    yplot_fit = int_rho_fit(xplot, *fit_coeffs)

    print(f"|| int_0^y rho dy - fit ||_oo = {la.norm(yplot_real - yplot_fit, np.inf)}")

    plt.plot(xplot, yplot_real, label="real")
    plt.plot(xplot, yplot_fit, label="fit")
    plt.plot(xplot, xplot, label="y", linestyle="--")
    plt.xlabel("x")
    plt.ylabel("int_0^y rho dy")
    plt.grid()
    plt.legend()

    # now we need to solve for y in the equation
    # etaMax = fit(y), whose solution is yMax

    eta_max = x[-1]
    alpha, beta, a, b = fit_coeffs

    def func_to_solve(yval):
        return (
            int_rho_fit(yval, alpha, beta, a, b) * param.deltaStar / deltaStar_inflow
            - eta_max
        )

    # initial guess: last y value (physical guess)
    y_max = fsolve(func_to_solve, y_physical_inflow[-1])[0]
    print(f"y corresponding to eta_max = {eta_max} is y_max = {y_max}")

    # correct the fit_coeffs to include the pre-factor deltaStar / deltaStar_inflow and add eta_max and y_max
    tmp = param.deltaStar / deltaStar_inflow
    fit_coeffs = np.array([tmp, tmp * alpha, beta, tmp * a, b, eta_max, y_max])

    return fit_coeffs


def printCoeffs_incNS(coeffs: list, param: Parameters, limits: list) -> None:
    # for rational function
    labels = ["u", "v"]
    for coeff, limit, label in zip(coeffs, limits, labels):
        for i, c in enumerate(coeff):
            if i < param.p:
                print(f"<p> a{i + 1}_{label} = {c} </p>")
            else:
                print(f"<p> b{i - param.p + 1}_{label} = {c} </p>")
        print(f"<p> b{param.p}_{label} = {coeff[param.p - 1] / limit} </p>")


def printCoeffs_comNS(coeffs: list, param: Parameters, limits: list) -> None:
    # for rational function
    labels = ["rho", "rhou", "rhov", "E1", "E2"]
    labels_y2eta = ["a0", "a1", "b1", "a2", "b2"]
    coeffs_y2eta = coeffs.pop()
    yMax = coeffs_y2eta[-1]
    etaMax = coeffs_y2eta[-2]
    coeffs_y2eta = coeffs_y2eta[:-2]

    print("")  # extra space
    print(f"""
    <p> C               = {param.deltaStar} </p> 
    <p> Ma              = {param.mainf} </p> <!-- REMEBER CHANGING BLASIUS PROFILE IF CHANGING THIS !-->
    <p> Pr              = {param.pr} </p> <!-- around 0.71 and 0.72 for air !-->
    <p> Gamma           = {param.gamma} </p>
    <p> Tref            = {param.tref_inf} </p> <!-- [K], reference temperature for Sutherland's law !-->
    <p> etaMax          = {etaMax} </p>
    <p> yMax            = {yMax} </p>
    """)

    print("")  # extra space

    for label, c in zip(labels_y2eta, coeffs_y2eta):
        print(f"      <p> {label}_y2eta = {c} </p>")
    for field_coeffs, limit, label in zip(coeffs, limits, labels):
        for i, c in enumerate(field_coeffs):
            if i < param.p:
                if label in ["rho", "rhou", "E1"]:
                    print(f"      <p> a{i + 1}_{label} = {c} </p>")
                else:
                    print(f"      <p> a{i + 2}_{label} = {c} </p>")
            else:
                print(f"      <p> b{i - param.p + 1}_{label} = {c} </p>")
        if label in ["rho", "rhou", "E1"]:
            print(
                f"      <p> b{param.p}_{label} = {field_coeffs[param.p - 1] / limit} </p>"
            )
        else:
            print(
                f"      <p> b{param.p + 1}_{label} = {field_coeffs[param.p - 1] / limit} </p>"
            )

    print("""    </PARAMETERS>

    <VARIABLES>
      <V ID="0"> rho  </V>
      <V ID="1"> rhou </V>
      <V ID="2"> rhov </V>
      <V ID="3"> E    </V>
    </VARIABLES>
     
    <BOUNDARYREGIONS>
      <B ID="0"> C[1] </B>    <!-- inlet !-->  
      <B ID="1"> C[2] </B>    <!-- outlet !-->
      <B ID="2"> C[3] </B>    <!-- top !-->
      <B ID="3"> C[4] </B>    <!-- bottom + left & right gap walls !-->
    </BOUNDARYREGIONS>    

    <BOUNDARYCONDITIONS>
      <REGION REF="0">  <!-- inlet !-->""")
    eta = (
        "a0_y2eta * y + a1_y2eta * tanh(b1_y2eta * y) + a2_y2eta * tanh(b2_y2eta * y)^2"
    )
    for label in labels[:3]:
        freestreamLabel = f"{label}Inf"
        onePlus = "1 + " if label == "rho" else ""
        add1 = 1 if label in ["rhov", "E2"] else 0
        print(f"""        <D VAR="{label}" VALUE="{freestreamLabel} * ( 
                   (
                     ({onePlus}{" + ".join([f"a{i + 1 + add1}_{label}*({eta})^{i + 1 + add1}" for i in range(param.p)])}) / 
                     (1 + {" + ".join([f"b{i + 1}_{label}*({eta})^{i + 1}" for i in range(param.p + add1)])})
                   ) * (y<yMax) + 
                   (
                     ({onePlus}{" + ".join([f"a{i + 1 + add1}_{label}*etaMax^{i + 1 + add1}" for i in range(param.p)])}) / 
                     (1 + {" + ".join([f"b{i + 1}_{label}*etaMax^{i + 1}" for i in range(param.p + add1)])})
                   ) * (y>=yMax))" />""")

    # energy equation
    print(f"""        <D VAR="E" VALUE="E0Inf + E1Inf * (
                   (
                     ({" + ".join([f"a{i + 1}_E1*({eta})^{i + 1}" for i in range(param.p)])}) / 
                     (1 + {" + ".join([f"b{i + 1}_E1*({eta})^{i + 1}" for i in range(param.p)])})
                   ) * (y<yMax) + 
                   (
                     ({" + ".join([f"a{i + 1}_E1*etaMax^{i + 1}" for i in range(param.p)])}) / 
                     (1 + {" + ".join([f"b{i + 1}_E1*etaMax^{i + 1}" for i in range(param.p)])})
                   ) * (y>=yMax)) + E2Inf * (
                   (
                     ({" + ".join([f"a{i + 2}_E2*({eta})^{i + 2}" for i in range(param.p)])}) / 
                     (1 + {" + ".join([f"b{i + 1}_E2*({eta})^{i + 1}" for i in range(param.p + 1)])})
                   ) * (y<yMax) + 
                   (
                     ({" + ".join([f"a{i + 2}_E2*etaMax^{i + 2}" for i in range(param.p)])}) / 
                     (1 + {" + ".join([f"b{i + 1}_E2*etaMax^{i + 1}" for i in range(param.p + 1)])})
                   ) * (y>=yMax))" />""")

    print("""      </REGION>""")
