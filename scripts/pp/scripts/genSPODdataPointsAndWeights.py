import numpy as np
import matplotlib.pyplot as plt
import os
from pp.fileManagement import extract_width_depth
from pp.inputargs import parseArgs
from scipy.spatial import Voronoi, voronoi_plot_2d

def generate_points_outside_gap(xmin: float, xmax: float, ymax: float, num_points_x: int, num_points_y: int) -> np.ndarray:
    x = np.linspace(xmin, xmax, num_points_x)
    y = np.linspace(0, ymax, num_points_y)
    X, Y = np.meshgrid(x, y)
    points = np.vstack([X.ravel(), Y.ravel()]).T

    return points


def generate_points_inside_gap(depth: float, width: float, num_points_x: int, num_points_y: int) -> np.ndarray:
    x = np.linspace(0, width, num_points_x)
    y = np.linspace(-depth, 0, num_points_y) [:-1]  # exclude the last point to avoid duplication at y=0
    X, Y = np.meshgrid(x, y)
    points = np.vstack([X.ravel(), Y.ravel()]).T
    return points

def print_points(points: np.ndarray, folder: str) -> None:
    filename = "points.txt"
    filename = os.path.join(folder, filename)
    with open(filename, 'w') as f:
        for point in points:
            f.write(f"{point[0]} {point[1]} 0\n")
    print(f"{len(points)} points written to {filename}")

def print_weights(weights: np.ndarray, folder: str) -> None:
    filename = "weights.txt"
    filename = os.path.join(folder, filename)
    with open(filename, 'w') as f:
        for weight in weights:
            f.write(f"{weight}\n")
    print(f"{len(weights)} weights written to {filename}")

def plot_points(points: np.ndarray) -> None:
    # rainbow colormap based on the index of the point
    # also plot a small label with the index at each point
    vor = Voronoi(points)
    # fig = voronoi_plot_2d(vor)


    scatter = plt.scatter(points[:, 0], points[:, 1],
                          c=np.arange(len(points)), cmap='gist_rainbow', s=5)

    # Add index labels next to each point
    freq = 10
    for i, (x, y) in enumerate(points):
        if i % freq == 0:
            plt.text(x, y, str(i), fontsize=6, color='black')

    plt.colorbar(scatter, label='Point Index')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Generated Points')
    plt.grid(True)
    plt.show()

def compute_weights(points_outside, points_inside):
    """
    Compute rectangular weights for points inside and outside the gap.
    The logic matches Victor's boundary definitions.
    """
    # --- Separate and sort grids ---
    # OUTSIDE region
    x_out = np.unique(points_outside[:, 0])
    y_out = np.unique(points_outside[:, 1])
    dx_out = np.gradient(x_out)
    dy_out = np.gradient(y_out)

    DX_out, DY_out = np.meshgrid(dx_out, dy_out)
    W_out = DX_out * DY_out

    # Apply boundary halves
    W_out[0, :]    /= 2   # bottom boundary (y=0)
    W_out[-1, :]   /= 2   # top boundary (y=h)
    W_out[:, 0]    /= 2   # left boundary (x=-Li)
    W_out[:, -1]   /= 2   # right boundary (x=Lo)

    # INSIDE region
    x_in = np.unique(points_inside[:, 0])
    y_in = np.unique(points_inside[:, 1])
    dx_in = np.gradient(x_in)
    dy_in = np.gradient(y_in)

    DX_in, DY_in = np.meshgrid(dx_in, dy_in)
    W_in = DX_in * DY_in

    # Apply boundary halves
    W_in[0, :]    /= 2    # bottom boundary (y=-d)
    W_in[:, 0]    /= 2    # left boundary (x=0)
    W_in[:, -1]   /= 2    # right boundary (x=w)
    W_in[-1, :]  *= 1.5  # top boundary (y=0)

    # Flatten weights back to match point order
    w_out = W_out.ravel()  # same order as points_outside
    w_in  = W_in.ravel()   # same order as points_inside

    return np.concatenate([w_out, w_in])


def main():
    
    # get path working directory
    args = parseArgs()

    folder = args.folders[0] # we do one by one

    d, w = extract_width_depth(folder)

    xmin = -10
    xmax = w + 4 * w
    ymax = 10

    dy = 0.2
    aspect_ratio = (xmax - xmin) / (ymax + d)
    dx = dy * aspect_ratio


    num_points_x_outside = int((xmax - xmin) / dx) + 1
    num_points_y_outside = int(ymax / dy) + 1
    points_outside = generate_points_outside_gap(xmin, xmax, ymax, num_points_x_outside, num_points_y_outside)

    num_points_x_inside = int(w / dx) + 1
    num_points_y_inside = int(d / dy) + 1
    points_inside = generate_points_inside_gap(d, w, num_points_x_inside, num_points_y_inside)

    all_points = np.vstack([points_outside, points_inside])
    weights = compute_weights(points_outside, points_inside)

    print("Those two quantities should be equal:")
    print(f"Total area from weights: {np.sum(weights)}")
    print(f"Geometrical area: { (xmax - xmin) * ymax + (w * d)}")
    
    # point_i = 480
    # print(f"Point {point_i}: {all_points[point_i]}")

    print("dx:", dx)
    print("dy:", dy)
    print_points(all_points, folder)
    print_weights(weights, folder)

    plot_points(all_points)



if __name__ == "__main__":
    main()
