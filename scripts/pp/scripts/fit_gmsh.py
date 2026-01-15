import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def main():
    """
    Main function to fit data points to a sigmoid-like function and plot the results.
    It is used in .geo files from gmsh to use a autmatied script depending on the width (x) of the gap.
    The yleft and yright values are the geometric rates in the pooints in the left part and right part of the gap, respectively.
    """
    # Data points
    x = np.array([5, 10, 16, 30, 50, 75, 100, 150])
    yleft = np.array([0.6, 0.76, 0.84, 0.89, 0.935, 0.96, 0.97, 0.98])
    yright = np.array([0.6, 0.8, 0.87, 0.93, 0.95, 0.97, 0.98, 0.985])

    # Define a function to fit (sigmoid-like function)
    def sigmoid(x, a, b, c, d):
        return c - a / (x - b) ** d

    # Fit both y1 and y2 using curve fitting
    params_yleft, _ = curve_fit(sigmoid, x, yleft, maxfev=10000)
    params_yright, _ = curve_fit(sigmoid, x, yright, maxfev=10000)

    # Print the parameters for both fits
    print(f"p_in_v_left = {params_yleft[2]} - {params_yleft[0]}/(W - {params_yleft[1]})^{params_yleft[3]}; // densitiy concentration of elements inside the gap, vertically, near the left wall.")
    print(f"p_in_v_right = {params_yright[2]} - {params_yright[0]}/(W - {params_yright[1]})^{params_yright[3]}; // densitiy concentration of elements inside the gap, vertically, near the right wall.")

    # Generate y-values from the fitted functions
    x_fit = np.linspace(min(x), max(x), 100)
    y1_fit = sigmoid(x_fit, *params_yleft)
    y2_fit = sigmoid(x_fit, *params_yright)

    # Plot the original data and the fitted functions
    plt.figure(figsize=(10, 6))
    plt.scatter(x, yleft, color="blue", label="y1 data")
    plt.scatter(x, yright, color="green", label="y2 data")
    plt.plot(x_fit, y1_fit, color="blue", linestyle="--", label="y1 fit")
    plt.plot(x_fit, y2_fit, color="green", linestyle="--", label="y2 fit")

    plt.xlabel("w")
    plt.ylabel("Geometric Rate")
    plt.title("Data and Sigmoid Fit")
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()
