import numpy as np
import matplotlib.pyplot as plt

# --- fixed geometry (inches) -------------------------------------------------
FRAME_W = 6.0  # width of the actual plot frame (the axes box)
FRAME_H = FRAME_W * 0.75  # height of the frame
PAD_X = 3.0  # space added to the LEFT and to the RIGHT of the frame
PAD_B = 0.8  # space below the frame (xlabel / xticks)
PAD_T = 0.4  # space above the frame (title)


def quick_plot():
    fig_w = FRAME_W + 2 * PAD_X
    fig_h = FRAME_H + PAD_B + PAD_T
    fig = plt.figure(figsize=(fig_w, fig_h))

    # place the axes by hand, in figure fractions: [left, bottom, width, height].
    # since the left and right pads are equal, the frame is exactly centered,
    # independently of the ylabel/tick widths or of the legend on the right.
    ax = fig.add_axes([PAD_X / fig_w, PAD_B / fig_h, FRAME_W / fig_w, FRAME_H / fig_h])

    x = [46.21, 46.32, 46.46, 46.88, 47.1, 47.5]
    y = [1.77, 1.78, 1.79, 1.8, 1.81, 1.82]

    ax.plot(x, y, marker="o", linestyle="-", color="b")
    x_fit = np.linspace(44.5, 48, 100)

    # get the regression line
    coeffs = np.polyfit(x, y, 1)
    y_fit = np.polyval(coeffs, x_fit)
    ax.plot(x_fit, y_fit, linestyle="--", color="g", label="Regression Line")
    print(f"Regression line equation: y = {coeffs[0]:.4f}x + {coeffs[1]:.4f}")

    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    ax.set_title("Quick Plot")
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.grid()

    fig.savefig("../../images/quick_plot.pdf", format="pdf")


if __name__ == "__main__":
    quick_plot()
