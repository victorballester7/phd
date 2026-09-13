import numpy as np
import matplotlib.pyplot as plt


def quick_plot():
    plt.figure(figsize=(10, 6))
    x = [46.21, 46.32, 46.46, 46.88, 47.1, 47.5]
    y = [1.77, 1.78, 1.79, 1.8, 1.81, 1.82]

    plt.plot(x, y, marker="o", linestyle="-", color="b")
    x_fit = np.linspace(44.5, 48, 100)


    # get the regression line
    coeffs = np.polyfit(x, y, 1)
    y_fit = np.polyval(coeffs, x_fit)
    plt.plot(x_fit, y_fit, linestyle="--", color="g", label="Regression Line")
    print(f"Regression line equation: y = {coeffs[0]:.4f}x + {coeffs[1]:.4f}")

    plt.title("Quick Plot")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    quick_plot()
