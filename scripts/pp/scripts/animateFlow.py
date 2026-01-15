import pyvista as pv
import os
import time
import re
import numpy as np
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial

class Bounds:
    def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax


def get_time_from_xml(xml_path: str) -> float:
    """Extract time value from Info.xml file."""
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        
        # Find the Time element in Metadata
        time_element = root.find('.//Metadata/Time')
        if time_element is not None:
            return float(time_element.text)
        else:
            print(f"Warning: Time element not found in {xml_path}")
            return None
    except Exception as e:
        print(f"Error reading {xml_path}: {e}")
        return None


def load_data(file_path: str):
    """Load VTU files and extract corresponding times from Info.xml files."""
    
    pattern = re.compile(r"mesh_(\d+)\.vtu")
    values = []
    
    for f in os.listdir(file_path):
        match = pattern.match(f)
        if match:
            num = int(match.group(1))
            values.append(num)

    values = np.array(values)
    min_val = np.min(values)
    max_val = np.max(values)

    print(f"Loading files from {min_val} to {max_val}")
        
    files = []
    times = []

    for f in range(min_val, max_val + 1):
        # Get time from Info.xml
        xml_path = os.path.join(file_path, f"mesh_{f}.chk", "Info.xml")
        print(f"Getting time for file mesh_{f}.chk")
        
        sim_time = get_time_from_xml(xml_path)
        if sim_time is not None:
            times.append(sim_time)
        else:
            # If time not found, use sequential numbering as fallback
            times.append(float(f))
        
        # Load VTU file
        print(f"Loading file mesh_{f}.vtu")
        files.append(pv.read(os.path.join(file_path, f"mesh_{f}.vtu")))

    print("All files loaded!")
    return files, np.array(times)


def clip_mesh(meshes: list, bounds: Bounds):
    """Clip meshes to specified bounds."""
    clipped_meshes = []
    bounds_tuple = (
        bounds.xmin,
        bounds.xmax,
        bounds.ymin,
        bounds.ymax,
        bounds.zmin,
        bounds.zmax,
    )
    for i, mesh in enumerate(meshes):
        print(f"Clipping mesh {i}/{len(meshes)}")
        surf = mesh.extract_surface()
        cropped = surf.clip_box(bounds_tuple, invert=False)
        clipped_meshes.append(cropped)
    print("All meshes clipped!")
    return clipped_meshes


def create_gif_with_timing(clipped_meshes, times, output_path, clim=(-0.035, 0.035),
                           time_scale=0.5, min_frame_duration=0.1):
    """
    Create GIF with frame durations proportional to simulation time differences.
    
    Args:
        clipped_meshes: List of PyVista meshes
        times: Array of simulation times
        output_path: Path to save the GIF
        time_scale: Scale factor for frame durations (1.0 = real time)
        min_frame_duration: Minimum duration for any frame (seconds)
    """
    # Calculate time differences between frames
    time_diffs = np.diff(times)
    
    # Normalize time differences to frame durations
    # Scale them to reasonable values (e.g., 0.1 to 1.0 seconds per frame)
    if len(time_diffs) > 0 and np.max(time_diffs) > 0:
        frame_durations = time_diffs / np.max(time_diffs) * time_scale
        frame_durations = np.maximum(frame_durations, min_frame_duration)
    else:
        frame_durations = np.ones(len(clipped_meshes) - 1) * 0.2
    
    print(f"\nTime statistics:")
    print(f"  Simulation time range: {times[0]:.2f} to {times[-1]:.2f}")
    print(f"  Time step min/max: {np.min(time_diffs):.4f} / {np.max(time_diffs):.4f}")
    print(f"  Frame duration min/max: {np.min(frame_durations):.3f}s / {np.max(frame_durations):.3f}s")
    
    # Create plotter
    pl = pv.Plotter(off_screen=True)
    
    # Since PyVista's open_gif doesn't support variable frame rates,
    # we'll duplicate frames to approximate variable timing
    gif_fps = 10  # Fixed FPS for the GIF
    frame_dt = 1.0 / gif_fps

    # check if output path contins "inc", otherwise set variable to "rhov"
    variable = "v" if "inc" in output_path else "rhov"
    
    pl.open_gif(output_path, fps=gif_fps)
    
    for i, mesh in enumerate(clipped_meshes):
        pl.clear()
        pl.add_mesh(mesh, scalars=variable, cmap="viridis", clim=clim)
        pl.view_xy()
        pl.enable_parallel_projection()

        # add number of current frame as title
        pl.add_text(f"Frame {i+1}/{len(clipped_meshes)} | Time: {times[i]:.2f}", font_size=10)
        
        # Calculate how many frames to write for this mesh
        if i < len(frame_durations):
            n_frames = max(1, int(frame_durations[i] / frame_dt))
        else:
            n_frames = 1
        
        # Write multiple frames for longer durations
        for _ in range(n_frames):
            pl.write_frame()
        
        print(f"Frame {i+1}/{len(clipped_meshes)} (time={times[i]:.2f}, duration={frame_durations[i] if i < len(frame_durations) else 0:.3f}s, frames={n_frames})")
    
    pl.close()
    print(f"\nGIF saved to: {output_path}")


def main():
    d = 2.5
    w = 23
    codename = f"d{d}_w{w}"
    typeOfRun = "comGapRe1000Ma0.2"
    script_path = os.path.dirname(os.path.abspath(__file__))
    directory = os.path.expanduser(
        f"~/hosts/hpc/PhD/runs/{typeOfRun}/baseflow/dns/{codename}/"
    )

    # Load meshes and times
    meshes, times = load_data(directory)
    
    # Clip meshes
    bound = Bounds(0, w + 150, -d, 20, -1e-6, 1e-6)
    clipped_meshes = clip_mesh(meshes, bound)

    # Create GIF with timing proportional to simulation time
    output_path = os.path.join(script_path, f"../../../videos/{codename}_{typeOfRun}.gif")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    create_gif_with_timing(
        clipped_meshes, 
        times, 
        output_path, 
        clim=(-0.02, 0.02),
        time_scale=0.4,  # Adjust this to speed up/slow down the animation. The lower the value, the faster the animation.
        min_frame_duration=0.05  # Minimum time per frame
    )


if __name__ == "__main__":
    main()
