import numpy as np

class Parameters:
    def __init__(
        self,
        eta_int_max: float,
        n: int,
        p: int,
        eta_interpolation_max: float,
        x_inflow: float,
        incNS: bool,
        uinf: float,
        rhoinf: float,
        re_deltaStar: float,
        mainf: float,
        pr: float,
        tinf: float,
        tref_inf: float,
        twall_dimensionless: float,
        gamma: float,
        adiabatic: bool,
    ) -> None:
        self.eta_int_max = eta_int_max
        self.n = n
        self.p = p
        self.eta_interpolation_max = eta_interpolation_max
        self.x_inflow = x_inflow
        self.incNS = incNS
        self.uinf = uinf
        self.rhoinf = rhoinf
        self.re_deltaStar = re_deltaStar
        self.mainf = mainf
        self.pr = pr
        self.tinf = tinf
        self.tref_inf = tref_inf
        self.twall_dimensionless = twall_dimensionless
        self.gamma = gamma
        self.adiabatic = adiabatic
        self.dim_system = 3 if self.incNS else 5

        self.x_interp = np.array([])
        self.y_interp = np.array([])
    
        self.delta = 0
        self.deltaStar = 0
        self.theta = 0
        self.shapeFactor = 0

    def __repr__(self):
        return (
            "Boundary layer quantities (the formulas are valid for both incNS and comNS)\n"
            f"δ          = {self.delta} * x / sqrt(Re_x)\n"
            f"δ*         = {self.deltaStar} * x / sqrt(Re_x)\n"
            f"θ          = {self.theta} * x / sqrt(Re_x)\n"
            f"δ/δ*       = {self.delta / self.deltaStar}\n"
            f"H (= δ*/θ) = {self.shapeFactor}\n"
        )

    def computeBLParams(self, x: np.ndarray, y: np.ndarray) -> None:
        self.delta = self.computeDelta(x, y[1])
        self.deltaStar = self.computeDeltaStar(x, y, self.incNS)
        self.theta = self.computeTheta(x, y[1])
        self.shapeFactor = self.deltaStar / self.theta

    def setxy_interp(self, x: np.ndarray, y: np.ndarray) -> None:
        mask = x <= self.eta_interpolation_max
        self.x_interp = x[mask]
        self.y_interp = y[:, mask]
        

    # BL thickness
    @staticmethod
    def computeDelta(x: np.ndarray, u: np.ndarray) -> float:
        # thinkness boundary layer (as first index of x_plot such that y_plot[1] > 0.99)
        limit = 0.99
        aux = np.abs(u - limit)
        idx = np.where(aux == np.min(aux))[0]
        if u[idx] < limit:
            idx1 = idx
            idx2 = idx + 1
        else:
            idx1 = idx - 1
            idx2 = idx

        L = u[idx2] - u[idx1]  # which is > 0
        d1 = limit - u[idx1]
        d2 = L - d1
        delta = ((d1 * x[idx2] + d2 * x[idx1]) / L)[0]
        return delta


    # displacement thickness
    @staticmethod
    def computeDeltaStar(x: np.ndarray, y: np.ndarray, incNS: bool) -> float:
        # to compute deltaStar, we need to integrate [1 - rho * u / (rhoInf * uInf)] from 0 to infinity (actually till etamax, as approximated)
        # in the incNS setting this means integrating (1-y[1]) from 0 to infinity (actually till etamax, as approximated)
        # in the comNS setting this means integrating (1-y[1]/y[3]) from 0 to infinity (actually till etamax, as approximated)

        # deltaStar = quad(lambda x: 1 - np.polyval(coeffs_fprime, x), 0, eta_max)[0]
        # print("deltaStar = ", deltaStar)
        # use numpy to compute the integral (they are basically the same, around 1e-7 difference for N = 10000)
        # deltaStar = np.trapezoid(1 - fprime, x)
        # print("deltaStar = ", deltaStar)

        # exact value is (eta - f(eta))|_0^inf = lim_{eta->inf} eta - f(eta) + f(0) = lim_{eta->inf} eta - f(eta)
        if incNS:
            f = y[0]
            deltaStar = x[-1] - f[-1]
        else:
            f = y[0]
            g = y[3]
            deltaStar = np.trapezoid(g , x) - f[-1]
        return deltaStar


    # momentum thickness
    @staticmethod
    def computeTheta(x: np.ndarray, u: np.ndarray) -> float:
        # to compute theta, we need to integrate [rho * u / (rhoInf * uInf) * (1 - u / uInf)] from 0 to infinity (actually till etamax, as approximated)
        # in th incNS setting this means integrating y[1](1-y[1]) from 0 to infinity (actually till etamax, as approximated)
        # in the comNS setting this means integrating y[1] / y[3] * (1 - y[1]) from 0 to infinity (actually till etamax, as approximated)
        theta = float(np.trapezoid(u * (1 - u), x))
        return theta


    def getRe_x(self, re_deltaStar_x: float) -> float:
        return (re_deltaStar_x / self.deltaStar) ** 2


    def y2eta_incNS(self, y: np.ndarray, deltaStar_x: float) -> np.ndarray:
        """
        Convert y to eta using the given parameters.
        """
        return y * self.deltaStar / deltaStar_x
