import argparse
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import numpy as np
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
from matplotlib.path import Path
from typing import Tuple
from pp.colors import colors

def read_vtu(file_path: str) -> Tuple[vtkUnstructuredGrid, float, float]:
    """Reads a .vtu file and returns the VTK unstructured grid object."""
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_path)
    reader.Update()

    # exctract the number X and Y of the file which ends with path/to/file/mesh_d4_w16.35_83.vtu
    file_path = file_path.split("/")[-1]
    d = float(file_path.split("_")[1][1:])
    w = float(file_path.split("_")[2][1:])

    return reader.GetOutput(), d, w

def colormap(field_name: str) -> str:
    """Returns the colormap to use for a given field."""
    if field_name in ["u", "v", "w"]:
        return "coolwarm"
    elif field_name in ["p"]:
        return "viridis"
    elif field_name in ["W_x", "W_y", "W_z"]:
        return "plasma"
    else:
        return "jet"

def getDomainBoundary(d: float, w: float) -> np.ndarray:
    """Returns the boundary of the domain."""
    
    lengthInflow = 100
    lengthOutflow = 1000
    lengthTop = 150
    
    boundary_points = np.array([
        [0, 0], [0, -d], [w, -d], [w, 0], [w + lengthOutflow, 0],
        [w + lengthOutflow, lengthTop], [-lengthInflow, lengthTop], [-lengthInflow, 0], [0, 0]
    ])

    return boundary_points

def setPlotLimits(x_coords: np.ndarray, y_coords: np.ndarray, d, w) -> Tuple[float, float, float, float]:
    """Sets the limits of the plot."""
    xmin, xmax = np.min(x_coords), np.max(x_coords)
    ymin, ymax = np.min(y_coords), np.max(y_coords)

    # xmin = 0
    # xmax = w
    # ymin = -d
    # ymax = 0

    return xmin, xmax, ymin, ymax


def plot_vtu(data: vtkUnstructuredGrid, d: float, w: float, field_name: str, save: str = "") -> None:
    """Plots a heatmap of a given field from the .vtu data while ensuring the domain boundary is respected."""
    point_data = data.GetPointData().GetArray(field_name)
    if point_data is None:
        print(f"Field '{field_name}' not found in the dataset.")
        return

    num_points = data.GetNumberOfPoints()
    coords = np.array([data.GetPoint(i) for i in range(num_points)])
    values = np.array([point_data.GetValue(i) for i in range(num_points)])

    x_coords, y_coords = coords[:, 0], coords[:, 1]

    # Define the boundary of the domain
    boundary_points = getDomainBoundary(d, w)
    boundary_path = Path(boundary_points)

    # Create a Delaunay triangulation
    print(colors.OKBLUE + "Creating Delaunay triangulation..." + colors.ENDC)
    triang = tri.Triangulation(x_coords, y_coords)

    # Mask triangles outside the boundary
    mask = np.array([not boundary_path.contains_point(np.mean(coords[triangle], axis=0))
                      for triangle in triang.triangles])
    triang.set_mask(mask)

    print(colors.OKBLUE + "Plotting..." + colors.ENDC)
    plt.figure()
    plt.tricontourf(triang, values, levels=100, cmap=colormap(field_name))

    cbar = plt.colorbar()
    cbar.ax.set_ylabel(field_name)

    xmin, xmax, ymin, ymax = setPlotLimits(x_coords, y_coords, d, w)
    # plt.plot(boundary_points[:, 0], boundary_points[:, 1], 'k-', linewidth=1.5)  # Draw the boundary
    # aspect_ratio set to 'equal' to ensure the domain is not distorted
    plt.gca().set_aspect('equal', adjustable='box')

    plt.axis("off")  # Remove figure framing
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    
    if save:
        plt.savefig(save, bbox_inches='tight')
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="Read and process a .vtu file.")
    parser.add_argument("file_path", type=str, help="Path to the .vtu file.")
    parser.add_argument("--variable", default="u", type=str, help="Variable to plot (default: u).")
    parser.add_argument("--save", default="", type=str, help="Path to save the plot (default: show only).")
    args = parser.parse_args()

    print(colors.OKBLUE + "Reading .vtu file..." + colors.ENDC)
    data, d, w = read_vtu(args.file_path)
    plot_vtu(data, d, w, args.variable, save=args.save)

if __name__ == "__main__":
    main()
