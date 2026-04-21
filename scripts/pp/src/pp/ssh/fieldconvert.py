import os
import xml.etree.ElementTree as ET

from pp.colors import colors
from pp.ssh.utilities import run_cmd, get_remote_dir


def fld2datapts(
    dir_local: str,
    fld_file: str,
    output_file: str,
    mesh: str = "mesh.xml",
    host: str = "hpc",
    user: str = "vb824",
):
    """
    Interpolates data from a .fld file to points defined in a .pts file using FieldConvert on a remote HPC cluster.

    Parameters:
    - dir_local: Local directory containing the session and mesh files (e.g., /home/Desktop/...../incNS/baseflow/dns/d3_w24/)
    - dir_remote: Remote directory on the HPC cluster where the interpolation will be performed (equivalent to dir_local but on the HPC)
    - fld_file: Name of the .fld file to be interpolated (e.g., mesh_35.chk)
    - output_file: Base name for the output files without extension (typically "data/filename" to store it in the data folder)
    - mesh: Name of the mesh file (default is "mesh.xml")
    - host: Hostname of the HPC cluster (default is "hpc")
    - user: Username for SSH access to the HPC cluster (default is "vb824")
    """
    print(colors.OKBLUE + "Starting interpolation process..." + colors.ENDC)

    dir_remote = get_remote_dir(dir_local)

    # get the parent directory of {DIR_REMOTE}/{output_file}.pts
    dir_remote_parent_output_file = os.path.dirname(os.path.join(dir_remote, output_file + ".pts"))

    # Copy pts file
    run_cmd(f"rsync -av {dir_local}/{output_file}.pts {user}@{host}:{dir_remote_parent_output_file}")

    # Run FieldConvert remotely
    remote_cmd = f"""
    source /etc/profile
    source ~/.bashrc  # Ensure modules are available
    cd {dir_remote}
    rm -rf {dir_remote_parent_output_file}/*.dat
    FieldConvert -m interppoints:fromxml={mesh}:fromfld={fld_file}:topts={output_file}.pts {output_file}.dat
    """

    run_cmd(f'ssh {user}@{host} "{remote_cmd}"')

    # Copy result back
    run_cmd(f"scp {user}@{host}:{dir_remote}/{output_file}.dat {dir_local}")

    print(
        colors.OKGREEN
        + "Interpolation process completed and data copied back to local machine."
        + colors.ENDC
    )


def combineAvg(
    dir_local: str,
    chkfile_ending_pattern: str,
    mintime: float,
    combine_avg: bool = True,
    mesh: str = "mesh.xml",
    host: str = "hpc",
    user: str = "vb824",
):
    """
    Combine averaged stress fields removing transient part.

    Uses mounted HPC directory for file inspection, only SSH for FieldConvert.

    Parameters:
    - dir_local:  Local directory containing the session and mesh files (e.g., /home/Desktop/...../incNS/baseflow/dns/d3_w24/)
    - dir_remote: Remote directory on the HPC cluster where the interpolation will be performed (equivalent to dir_local but on the HPC)
    - chkfile: suffix of stress folders
    - mesh: mesh file name
    - combine_avg: whether to perform combineAvg
    - mintime: minimum FinalTime threshold
    """

    final_avg_file = "mesh_avg.fld"

    print(colors.OKBLUE + "Scanning stress folders..." + colors.ENDC)

    dir_remote = get_remote_dir(dir_local)

    # Convert remote path → mounted path
    dir_mounted = dir_remote.replace(
        f"/rds/general/user/{user}/home/Desktop",
        f"/home/victor/hosts/{host}"
    )

    if not os.path.isdir(dir_mounted):
        raise FileNotFoundError(colors.FAIL + f"Mounted directory not found: {dir_mounted}" + colors.ENDC)

    # -----------------------
    # Find stress folders
    # -----------------------
    stress_folders = []

    for d in os.listdir(dir_mounted):
        if d.endswith(chkfile_ending_pattern):
            fld_path = os.path.join(dir_mounted, d, "P0000000.fld")
            if os.path.exists(fld_path):
                mtime = os.path.getmtime(fld_path)
                stress_folders.append((mtime, d))

    if not stress_folders:
        raise RuntimeError(colors.FAIL + f"No stress folders found with suffix '{chkfile_ending_pattern}' in {dir_mounted}" + colors.ENDC)

    # Sort by time
    stress_folders.sort()
    stress_folders = [d for _, d in stress_folders]

    stressFINAL = stress_folders[-1]
    print(colors.OKGREEN + f"Final folder: {stressFINAL}" + colors.ENDC)

    # -----------------------
    # If no combineAvg
    # -----------------------
    if not combine_avg:
        print(colors.WARNING + "Skipping combineAvg → copying final folder" + colors.ENDC)

        src = os.path.join(dir_mounted, stressFINAL)
        dst = os.path.join(dir_mounted, final_avg_file)

        run_cmd(f"rm -rf {dst}; cp -r {src} {dst}")
        return

    # -----------------------
    # Find start folder
    # -----------------------
    stressBEG = None

    for folder in stress_folders:
        xmlfile = os.path.join(dir_mounted, folder, "Info.xml")

        if not os.path.isfile(xmlfile):
            continue

        try:
            tree = ET.parse(xmlfile)
            root = tree.getroot()
            final_time = float(root.findtext(".//FinalTime", default="0"))
        except Exception:

            final_time = 0

        if final_time >= mintime:
            stressBEG = folder
            break

    if stressBEG is None:
        raise RuntimeError(colors.FAIL + f"No stress folder found with FinalTime >= {mintime}" + colors.ENDC)

    print(colors.OKGREEN + f"Start folder: {stressBEG}" + colors.ENDC)

    # -----------------------
    # Patch NumberOfFieldDumps
    # -----------------------
    xmlfile = os.path.join(dir_mounted, stressBEG, "Info.xml")

    try:
        tree = ET.parse(xmlfile)
        root = tree.getroot()

        node = root.find(".//NumberOfFieldDumps")
        if node is not None and node.text.isdigit():
            val = int(node.text)
            node.text = str(-val)
            tree.write(xmlfile)
            print(colors.WARNING + f"Patched NumberOfFieldDumps → {-val}" + colors.ENDC)

    except Exception as e:
        print(colors.WARNING + f"Failed to patch XML: {e}" + colors.ENDC)

    # -----------------------
    # Run FieldConvert remotely
    # -----------------------
    print(colors.OKBLUE + "Running FieldConvert combineAvg..." + colors.ENDC)

    remote_cmd = f"""
    source /etc/profile
    source ~/.bashrc
    cd {dir_remote}
    rm -rf mesh_avg.fld
    FieldConvert -m combineAvg:fromfld={stressBEG} {mesh} {stressFINAL} {final_avg_file}
    """

    run_cmd(f'ssh {user}@{host} "{remote_cmd}"')

    print(colors.OKGREEN + "combineAvg completed successfully." + colors.ENDC)

