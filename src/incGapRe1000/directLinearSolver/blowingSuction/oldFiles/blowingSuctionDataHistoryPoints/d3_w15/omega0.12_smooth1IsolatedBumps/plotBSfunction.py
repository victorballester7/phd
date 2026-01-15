import numpy as np
import matplotlib.pyplot as plt

def mollifier(x: np.ndarray, w: float) -> np.ndarray:
    """
    Mollifier function for smoothing.
    Args:
        x: Input array.
        w: Width of the mollifier.
    Returns:
        Smoothed array.
    """
    return np.exp(-w/ (1 - x**2)) * (np.abs(x) < 1)

def main():
    # Example data
    omega = 0.12
    w = 0.025
    T = 2 * np.pi / omega
    n = 20
    m = 1
    totPeriods = m + n
    numTotPeriods = 5
    times = np.linspace(0, numTotPeriods*totPeriods * T,10000)

    # xi(t) = t mod (totPeriods * T)
    xi_t = times % (totPeriods * T)
    phi_t = mollifier(2 * xi_t / (m * T) - 1, w)
    f_t = np.sin(omega * xi_t) * phi_t
    
    plt.figure(figsize=(10, 6))
    plt.plot(times, f_t, label='f(t)', color='blue')

    plt.grid()
    plt.xlabel('Time')
    plt.ylabel('f(t)')

    plt.show()



if __name__ == "__main__":
    main() 
