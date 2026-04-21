#!/usr/bin/env python3
"""Lambda-2 Criterion Visualization for Vortex Identification.

This script provides publication-ready visualization of vortex structures
using the Lambda-2 criterion (Jeong & Hussain, 1995). The Lambda-2 criterion
identifies vortex cores as connected regions where λ₂ < 0, where λ₂ is the
second eigenvalue of S² + Ω² (S: strain rate tensor, Ω: vorticity tensor).

Features:
    - 3D isosurface rendering of vortex structures
    - 2D slice visualization at specified locations
    - Publication-quality figure export
    - Interactive and batch processing modes

Usage:
    python L2Criterion.py <path_to_vtu_file> [options]

Example:
    python L2Criterion.py mesh_l2.vtu --threshold -0.01 --save figures/

References:
    Jeong, J., & Hussain, F. (1995). On the identification of a vortex.
    Journal of Fluid Mechanics, 285, 69-94.

Author: Victor
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib.colors import Normalize

# Try to load custom style if available
STYLE_PATH = Path(__file__).parent.parent / "plots" / "style" / "tsfp.mplstyle"
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))
else:
    plt.rcParams.update({
        "font.size": 12,
        "axes.labelsize": 14,
        "axes.titlesize": 14,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })


def load_lambda2_mesh(filepath: str, field_name: str = "L2") -> pv.DataSet:
    """Load mesh with Lambda-2 field data from VTU/VTK file.

    Args:
        filepath: Path to the VTU or VTK mesh file.
        field_name: Name of the Lambda-2 field in the mesh (default: "L2").

    Returns:
        PyVista mesh object with Lambda-2 data.

    Raises:
        FileNotFoundError: If the mesh file does not exist.
        KeyError: If the Lambda-2 field is not found in the mesh.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Mesh file not found: {filepath}")

    mesh = pv.read(filepath)

    if field_name not in mesh.array_names:
        available = ", ".join(mesh.array_names)
        raise KeyError(
            f"Field '{field_name}' not found. Available fields: {available}"
        )

    mesh.set_active_scalars(field_name)
    return mesh


def compute_lambda2_statistics(mesh: pv.DataSet, field_name: str = "L2") -> dict:
    """Compute statistics of the Lambda-2 field.

    Args:
        mesh: PyVista mesh with Lambda-2 data.
        field_name: Name of the Lambda-2 field.

    Returns:
        Dictionary with min, max, mean, std, and vortex volume fraction.
    """
    l2_values = mesh[field_name]

    stats = {
        "min": float(np.min(l2_values)),
        "max": float(np.max(l2_values)),
        "mean": float(np.mean(l2_values)),
        "std": float(np.std(l2_values)),
        "vortex_fraction": float(np.sum(l2_values < 0) / len(l2_values)),
    }

    return stats


def plot_lambda2_isosurface(
    mesh: pv.DataSet,
    threshold: float = -0.01,
    field_name: str = "L2",
    cmap: str = "coolwarm",
    opacity: float = 0.8,
    show_edges: bool = False,
    background: str = "white",
    camera_position: Optional[str] = "iso",
    window_size: Tuple[int, int] = (1920, 1080),
    save_path: Optional[str] = None,
    show: bool = True,
) -> pv.Plotter:
    """Render 3D isosurface of vortex structures using Lambda-2 criterion.

    Vortex cores are identified as regions where λ₂ < threshold.

    Args:
        mesh: PyVista mesh with Lambda-2 data.
        threshold: Lambda-2 threshold for isosurface (default: -0.01).
        field_name: Name of the Lambda-2 field.
        cmap: Colormap for the isosurface.
        opacity: Surface opacity (0-1).
        show_edges: Whether to show mesh edges.
        background: Background color.
        camera_position: Camera view ('iso', 'xy', 'xz', 'yz', or tuple).
        window_size: Render window size in pixels.
        save_path: Path to save the rendered image.
        show: Whether to display the interactive plot.

    Returns:
        PyVista Plotter object.
    """
    plotter = pv.Plotter(window_size=window_size, off_screen=not show)
    plotter.set_background(background)

    # Extract isosurface at threshold
    contours = mesh.contour(
        isosurfaces=[threshold],
        scalars=field_name,
        compute_normals=True,
        compute_scalars=True,
    )

    if contours.n_points > 0:
        plotter.add_mesh(
            contours,
            scalars=field_name,
            cmap=cmap,
            opacity=opacity,
            show_edges=show_edges,
            smooth_shading=True,
            scalar_bar_args={
                "title": r"$\lambda_2$",
                "title_font_size": 16,
                "label_font_size": 14,
                "shadow": True,
                "fmt": "%.2e",
                "position_x": 0.85,
                "position_y": 0.3,
                "height": 0.4,
            },
        )
    else:
        print(f"Warning: No isosurface found at threshold {threshold}")

    # Add domain outline
    plotter.add_mesh(mesh.outline(), color="black", line_width=2)

    # Set camera position
    if camera_position == "iso":
        plotter.camera_position = "iso"
        plotter.camera.azimuth = 30
        plotter.camera.elevation = 20
    elif camera_position in ("xy", "xz", "yz"):
        plotter.view_vector = {
            "xy": (0, 0, -1),
            "xz": (0, -1, 0),
            "yz": (-1, 0, 0),
        }[camera_position]
    elif isinstance(camera_position, (list, tuple)):
        plotter.camera_position = camera_position

    plotter.add_axes(
        xlabel="x",
        ylabel="y",
        zlabel="z",
        line_width=2,
        labels_off=False,
    )

    if save_path:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        plotter.screenshot(save_path, transparent_background=False)
        print(f"Saved 3D visualization to: {save_path}")

    if show:
        plotter.show()

    return plotter


def plot_lambda2_slices(
    mesh: pv.DataSet,
    slice_positions: Optional[dict] = None,
    field_name: str = "L2",
    cmap: str = "RdBu_r",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    symmetric: bool = True,
    figsize: Tuple[float, float] = (12, 8),
    save_path: Optional[str] = None,
    show: bool = True,
) -> plt.Figure:
    """Create 2D slice plots of the Lambda-2 field.

    Args:
        mesh: PyVista mesh with Lambda-2 data.
        slice_positions: Dict with 'x', 'y', 'z' keys specifying slice locations.
                         Default: slices at domain center.
        field_name: Name of the Lambda-2 field.
        cmap: Colormap for contour plots.
        vmin: Minimum value for colormap.
        vmax: Maximum value for colormap.
        symmetric: If True, use symmetric colormap around zero.
        figsize: Figure size in inches.
        save_path: Path to save the figure.
        show: Whether to display the plot.

    Returns:
        Matplotlib Figure object.
    """
    bounds = mesh.bounds
    center = mesh.center

    # Default slice positions at domain center
    if slice_positions is None:
        slice_positions = {
            "x": center[0],
            "y": center[1],
            "z": center[2],
        }

    l2_data = mesh[field_name]
    if symmetric:
        abs_max = max(abs(np.min(l2_data)), abs(np.max(l2_data)))
        vmin = vmin or -abs_max
        vmax = vmax or abs_max
    else:
        vmin = vmin or np.min(l2_data)
        vmax = vmax or np.max(l2_data)

    norm = Normalize(vmin=vmin, vmax=vmax)

    # Create slices
    slices = {}
    slice_configs = [
        ("x", "normal_x", slice_positions.get("x"), ("y", "z"), (1, 2)),
        ("y", "normal_y", slice_positions.get("y"), ("x", "z"), (0, 2)),
        ("z", "normal_z", slice_positions.get("z"), ("x", "y"), (0, 1)),
    ]

    for name, normal, pos, axes_labels, axes_idx in slice_configs:
        if pos is not None:
            origin = list(center)
            origin[axes_idx[0] if name == "x" else (0 if name != "x" else 1)] = pos
            if name == "x":
                origin[0] = pos
            elif name == "y":
                origin[1] = pos
            else:
                origin[2] = pos
            try:
                slices[name] = mesh.slice(normal=normal, origin=origin)
            except Exception:
                slices[name] = None

    # Create figure
    n_slices = sum(1 for s in slices.values() if s is not None and s.n_points > 0)
    if n_slices == 0:
        print("Warning: No valid slices could be extracted")
        return None

    fig, axes = plt.subplots(1, n_slices, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    ax_idx = 0
    for name, normal, pos, axes_labels, axes_idx in slice_configs:
        slice_mesh = slices.get(name)
        if slice_mesh is None or slice_mesh.n_points == 0:
            continue

        ax = axes[ax_idx]
        points = slice_mesh.points
        values = slice_mesh[field_name]

        # Project points to 2D
        x_coord = points[:, axes_idx[0]]
        y_coord = points[:, axes_idx[1]]

        # Triangulation for unstructured data
        from matplotlib.tri import Triangulation
        try:
            triang = Triangulation(x_coord, y_coord)
            tcf = ax.tricontourf(
                triang, values, levels=50, cmap=cmap, norm=norm, extend="both"
            )
            ax.tricontour(
                triang, values, levels=[0], colors="k", linewidths=1.5, linestyles="-"
            )
        except Exception:
            # Fallback to scatter plot
            scatter = ax.scatter(
                x_coord, y_coord, c=values, cmap=cmap, norm=norm, s=1
            )

        ax.set_xlabel(f"${axes_labels[0]}$")
        ax.set_ylabel(f"${axes_labels[1]}$")
        ax.set_title(f"${name} = {pos:.3f}$")
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.3)

        ax_idx += 1

    # Add colorbar
    cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=axes,
        orientation="vertical",
        fraction=0.02,
        pad=0.02,
    )
    cbar.set_label(r"$\lambda_2$")

    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Saved 2D slices to: {save_path}")

    if show:
        plt.show()

    return fig


def plot_lambda2_volume_rendering(
    mesh: pv.DataSet,
    field_name: str = "L2",
    opacity_unit_distance: float = 1.0,
    cmap: str = "coolwarm",
    window_size: Tuple[int, int] = (1920, 1080),
    save_path: Optional[str] = None,
    show: bool = True,
) -> pv.Plotter:
    """Create volume rendering of the Lambda-2 field.

    This provides a semi-transparent view of the entire field,
    highlighting regions of strong negative Lambda-2.

    Args:
        mesh: PyVista mesh with Lambda-2 data.
        field_name: Name of the Lambda-2 field.
        opacity_unit_distance: Controls opacity falloff with distance.
        cmap: Colormap for volume rendering.
        window_size: Render window size in pixels.
        save_path: Path to save the rendered image.
        show: Whether to display the interactive plot.

    Returns:
        PyVista Plotter object.
    """
    plotter = pv.Plotter(window_size=window_size, off_screen=not show)
    plotter.set_background("white")

    # Create custom opacity transfer function (emphasize negative values)
    l2_data = mesh[field_name]
    l2_min, l2_max = np.min(l2_data), np.max(l2_data)

    opacity = [1.0, 0.6, 0.2, 0.0, 0.0, 0.0, 0.0]  # High opacity for negative values

    plotter.add_volume(
        mesh,
        scalars=field_name,
        cmap=cmap,
        opacity=opacity,
        opacity_unit_distance=opacity_unit_distance,
        shade=True,
        scalar_bar_args={
            "title": r"$\lambda_2$",
            "title_font_size": 16,
            "label_font_size": 14,
            "fmt": "%.2e",
        },
    )

    plotter.add_mesh(mesh.outline(), color="black", line_width=2)
    plotter.camera_position = "iso"
    plotter.add_axes(xlabel="x", ylabel="y", zlabel="z")

    if save_path:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        plotter.screenshot(save_path, transparent_background=False)
        print(f"Saved volume rendering to: {save_path}")

    if show:
        plotter.show()

    return plotter


def create_publication_figure(
    mesh: pv.DataSet,
    threshold: float = -0.01,
    field_name: str = "L2",
    slice_positions: Optional[dict] = None,
    figsize: Tuple[float, float] = (16, 10),
    save_path: Optional[str] = None,
    show: bool = True,
) -> plt.Figure:
    """Create a comprehensive publication figure with multiple views.

    Generates a multi-panel figure with:
    - 3D isosurface view
    - Multiple 2D slice views
    - Statistics annotation

    Args:
        mesh: PyVista mesh with Lambda-2 data.
        threshold: Lambda-2 threshold for isosurface.
        field_name: Name of the Lambda-2 field.
        slice_positions: Dict specifying slice locations.
        figsize: Figure size in inches.
        save_path: Path to save the figure.
        show: Whether to display the plot.

    Returns:
        Matplotlib Figure object.
    """
    # Compute statistics
    stats = compute_lambda2_statistics(mesh, field_name)

    # Get domain info
    bounds = mesh.bounds
    center = mesh.center

    if slice_positions is None:
        slice_positions = {"x": center[0], "y": center[1], "z": center[2]}

    # Create figure with custom layout
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 1])

    # 3D isosurface (rendered off-screen, then embedded)
    ax_3d = fig.add_subplot(gs[0, :])

    # Render 3D view off-screen
    plotter = pv.Plotter(off_screen=True, window_size=(1600, 800))
    plotter.set_background("white")

    contours = mesh.contour(
        isosurfaces=[threshold],
        scalars=field_name,
        compute_normals=True,
    )

    if contours.n_points > 0:
        plotter.add_mesh(
            contours,
            scalars=field_name,
            cmap="coolwarm",
            opacity=0.85,
            smooth_shading=True,
            show_scalar_bar=False,
        )

    plotter.add_mesh(mesh.outline(), color="black", line_width=2)
    plotter.camera_position = "iso"
    plotter.camera.azimuth = 30
    plotter.camera.elevation = 20

    # Capture screenshot
    img = plotter.screenshot(return_img=True)
    plotter.close()

    ax_3d.imshow(img)
    ax_3d.axis("off")
    ax_3d.set_title(
        rf"Vortex structures ($\lambda_2 = {threshold:.2e}$)",
        fontsize=14,
        fontweight="bold",
    )

    # Add statistics text
    stats_text = (
        f"$\\lambda_{{2,\\min}} = {stats['min']:.2e}$\n"
        f"$\\lambda_{{2,\\max}} = {stats['max']:.2e}$\n"
        f"Vortex fraction: {stats['vortex_fraction']:.1%}"
    )
    ax_3d.text(
        0.02, 0.98, stats_text,
        transform=ax_3d.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    # 2D slices
    l2_data = mesh[field_name]
    abs_max = max(abs(np.min(l2_data)), abs(np.max(l2_data)))
    norm = Normalize(vmin=-abs_max, vmax=abs_max)
    cmap = "RdBu_r"

    slice_configs = [
        ("x", slice_positions.get("x"), ("y", "z"), (1, 2)),
        ("y", slice_positions.get("y"), ("x", "z"), (0, 2)),
        ("z", slice_positions.get("z"), ("x", "y"), (0, 1)),
    ]

    for idx, (name, pos, axes_labels, axes_idx) in enumerate(slice_configs):
        ax = fig.add_subplot(gs[1, idx])

        if pos is None:
            ax.text(0.5, 0.5, "No slice", ha="center", va="center")
            continue

        origin = list(center)
        origin[["x", "y", "z"].index(name)] = pos

        try:
            slice_mesh = mesh.slice(normal=f"normal_{name}", origin=origin)

            if slice_mesh.n_points > 0:
                points = slice_mesh.points
                values = slice_mesh[field_name]

                x_coord = points[:, axes_idx[0]]
                y_coord = points[:, axes_idx[1]]

                from matplotlib.tri import Triangulation
                triang = Triangulation(x_coord, y_coord)

                tcf = ax.tricontourf(
                    triang, values, levels=50, cmap=cmap, norm=norm, extend="both"
                )
                ax.tricontour(
                    triang, values, levels=[0], colors="k", linewidths=1.0
                )

        except Exception as e:
            ax.text(0.5, 0.5, f"Error: {e}", ha="center", va="center", fontsize=8)

        ax.set_xlabel(f"${axes_labels[0]}$")
        ax.set_ylabel(f"${axes_labels[1]}$")
        ax.set_title(f"${name} = {pos:.2f}$")
        ax.set_aspect("equal")

    # Add colorbar
    cbar_ax = fig.add_axes([0.92, 0.1, 0.015, 0.35])
    cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=cbar_ax,
    )
    cbar.set_label(r"$\lambda_2$")

    if save_path:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"Saved publication figure to: {save_path}")

    if show:
        plt.show()

    return fig


def main():
    """Main entry point for Lambda-2 visualization script."""
    parser = argparse.ArgumentParser(
        description="Visualize vortex structures using the Lambda-2 criterion.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s mesh.vtu
  %(prog)s mesh.vtu --threshold -0.001 --save output/
  %(prog)s mesh.vtu --mode slices --slice-x 5.0 --slice-z 0.5
  %(prog)s mesh.vtu --mode publication --save figures/vortices.png
        """,
    )

    parser.add_argument(
        "mesh_file",
        type=str,
        help="Path to the VTU/VTK mesh file with Lambda-2 data.",
    )
    parser.add_argument(
        "--field", "-f",
        type=str,
        default="L2",
        help="Name of the Lambda-2 field in the mesh (default: L2).",
    )
    parser.add_argument(
        "--threshold", "-t",
        type=float,
        default=-0.01,
        help="Lambda-2 threshold for isosurface (default: -0.01).",
    )
    parser.add_argument(
        "--mode", "-m",
        choices=["isosurface", "slices", "volume", "publication"],
        default="isosurface",
        help="Visualization mode (default: isosurface).",
    )
    parser.add_argument(
        "--slice-x",
        type=float,
        default=None,
        help="X-position for slice (uses center if not specified).",
    )
    parser.add_argument(
        "--slice-y",
        type=float,
        default=None,
        help="Y-position for slice (uses center if not specified).",
    )
    parser.add_argument(
        "--slice-z",
        type=float,
        default=None,
        help="Z-position for slice (uses center if not specified).",
    )
    parser.add_argument(
        "--cmap", "-c",
        type=str,
        default="coolwarm",
        help="Colormap for visualization (default: coolwarm).",
    )
    parser.add_argument(
        "--opacity",
        type=float,
        default=0.8,
        help="Opacity for 3D rendering (0-1, default: 0.8).",
    )
    parser.add_argument(
        "--save", "-s",
        type=str,
        default=None,
        help="Path to save the output figure/image.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Don't display the interactive plot (useful for batch processing).",
    )
    parser.add_argument(
        "--stats",
        action="store_true",
        help="Print Lambda-2 field statistics.",
    )

    args = parser.parse_args()

    # Load mesh
    print(f"Loading mesh from: {args.mesh_file}")
    try:
        mesh = load_lambda2_mesh(args.mesh_file, args.field)
    except (FileNotFoundError, KeyError) as e:
        print(f"Error: {e}")
        sys.exit(1)

    print(f"  Mesh points: {mesh.n_points:,}")
    print(f"  Mesh cells:  {mesh.n_cells:,}")
    print(f"  Bounds: {mesh.bounds}")

    # Print statistics if requested
    if args.stats:
        stats = compute_lambda2_statistics(mesh, args.field)
        print("\nLambda-2 Statistics:")
        print(f"  Min:            {stats['min']:.4e}")
        print(f"  Max:            {stats['max']:.4e}")
        print(f"  Mean:           {stats['mean']:.4e}")
        print(f"  Std:            {stats['std']:.4e}")
        print(f"  Vortex fraction: {stats['vortex_fraction']:.2%}")

    # Build slice positions
    slice_positions = {}
    if args.slice_x is not None:
        slice_positions["x"] = args.slice_x
    if args.slice_y is not None:
        slice_positions["y"] = args.slice_y
    if args.slice_z is not None:
        slice_positions["z"] = args.slice_z

    if not slice_positions:
        slice_positions = None  # Use defaults

    # Determine save path
    save_path = args.save
    if save_path and os.path.isdir(save_path):
        base_name = Path(args.mesh_file).stem
        ext = ".png"
        save_path = os.path.join(save_path, f"{base_name}_lambda2_{args.mode}{ext}")

    show = not args.no_show

    # Run visualization
    print(f"\nGenerating {args.mode} visualization...")

    if args.mode == "isosurface":
        plot_lambda2_isosurface(
            mesh,
            threshold=args.threshold,
            field_name=args.field,
            cmap=args.cmap,
            opacity=args.opacity,
            save_path=save_path,
            show=show,
        )

    elif args.mode == "slices":
        plot_lambda2_slices(
            mesh,
            slice_positions=slice_positions,
            field_name=args.field,
            cmap="RdBu_r",
            save_path=save_path,
            show=show,
        )

    elif args.mode == "volume":
        plot_lambda2_volume_rendering(
            mesh,
            field_name=args.field,
            cmap=args.cmap,
            save_path=save_path,
            show=show,
        )

    elif args.mode == "publication":
        create_publication_figure(
            mesh,
            threshold=args.threshold,
            field_name=args.field,
            slice_positions=slice_positions,
            save_path=save_path,
            show=show,
        )

    print("Done.")


if __name__ == "__main__":
    main()

