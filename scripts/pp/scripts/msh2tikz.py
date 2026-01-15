import meshio
import numpy as np
import os
from pathlib import Path

def readfilter_mesh(meshFile: Path, limits: tuple) -> tuple:
    xmin, xmax, ymin, ymax = limits

    # --- Read mesh ---
    mesh = meshio.read(meshFile)
    points = mesh.points

    # --- Filter points by x/y range ---
    in_box = (
        (points[:, 0] >= xmin) & (points[:, 0] <= xmax) &
        (points[:, 1] >= ymin) & (points[:, 1] <= ymax)
    )
    filtered_indices = np.where(in_box)[0]
    filtered_set = set(filtered_indices)

    # --- Get triangle cells ---
    triangles = None
    for cell_block in mesh.cells:
        if cell_block.type == "triangle":
            triangles = cell_block.data
            break

    if triangles is None:
        raise ValueError("No triangular cells found in the mesh.")

    # --- Generate TikZ draw commands ---
    draw_commands = set()  # use a set to avoid duplicate edges

    return points, triangles, filtered_set, draw_commands

def postprocess_tikz(tikzFile: Path, points: np.ndarray, triangles: np.ndarray, filtered_set: set, draw_commands: set) -> None:
    drawn_edges = set()  # to avoid duplicate edges
    for tri in triangles:
        if any(v in filtered_set for v in tri):
            for i, j in [(0, 1), (1, 2), (2, 0)]:
                vi, vj = tri[i], tri[j]
                edge = tuple(sorted((vi, vj)))  # ensure consistent order
                if edge in drawn_edges:
                    continue  # already added

                drawn_edges.add(edge)

                x1, y1 = points[vi][:2]
                x2, y2 = points[vj][:2]
                cmd = f"\\draw[thin,gray] ({x1},{y1}) -- ({x2},{y2});"
                draw_commands.add(cmd)

    # --- Print output ---
    with open(tikzFile, "w") as f:
        for cmd in sorted(draw_commands):
        # print in a file data/tikz_commands.txt
            f.write(cmd + "\n")


def main():
    """
    Main function to convert a mesh file (.msh from gmsh) to TikZ draw commands for LaTeX. The script generates an external .tex file that then will be included in the main LaTeX document, in order to reduce loading time (the file can be quite large).
    """
    pathCurrentScript = os.path.dirname(os.path.abspath(__file__))

    meshFile = "../../../src/incNSboeingGapRe1000/baseflow/dns/d4_w15/mesh_d4_w15.msh"
    tikzFile = "../../../latex/Images/data/edgesMesh.tex"

    meshFile = Path(os.path.join(pathCurrentScript, meshFile))
    tikzFile = Path(os.path.join(pathCurrentScript, tikzFile))

    # --- Parameters ---
    xmin, xmax = -15, 30
    ymin, ymax = 0.25, 11

    # --- Read and filter mesh ---
    points, triangles, filtered_set, draw_commands = readfilter_mesh(meshFile, (xmin, xmax, ymin, ymax))

    # --- Postprocess TikZ commands ---
    postprocess_tikz(tikzFile, points, triangles, filtered_set, draw_commands)
    
    
if __name__ == "__main__":
    main()
