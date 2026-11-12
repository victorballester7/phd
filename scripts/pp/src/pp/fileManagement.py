from typing import Dict
import gzip
import re
import os
import numpy as np
from pp.colors import colors
from typing import Tuple, List, Dict


def extract_depth_width(filenameFolder: str) -> Tuple[float, float]:
    """
    Extract width (YYY) and depth (XXX) from the filename of the form *dXXX_wYYY*
    Args:
        filenameFolder: Path to the folder containing the file.
    Returns:
        Width and depth values extracted from the filename. Defaults to 0 if not found.
    """

    match = re.search(r"d(\d+(?:\.\d+)?)_w(\d+(?:\.\d+)?)", filenameFolder)
    depth = 0
    width = 0
    if match:
        depth = float(match.group(1))
        width = float(match.group(2))
    if not match:
        match = re.search(r"d(\d+(?:\.\d+)?)", filenameFolder)
        if match:
            depth = float(match.group(1))
            width = 0
        else:
            print(
                colors.WARNING
                + f"Could not extract width from filename: {filenameFolder}. Returning default values : depth = {depth}, width = {width}"
                + colors.ENDC
            )
    return float(depth), float(width)

def extract_reynolds(path: str) -> float:
    match = re.search(r"Re(\d+(?:\.\d+)?)", path)
    if match is None:
        raise ValueError(f"Could not extract Reynolds number from path: {path}")
    return float(match.group(1))


def historyParts(folder: str, name: str) -> List[str]:
    """
    Files that hold the history/energy output ``name`` (e.g. "HistoryPoints.his") of a run directory, oldest first.

    Since 2026-11-10 pbspro.sh keeps every finished run in ``name.gz`` (one header, then the data of all runs) and only
    the run in progress in plain ``name``; ``name.old`` is the previous scheme (still read if present). ``folder`` may
    be the directory or any of these files.
    """
    if os.path.isdir(folder):
        base = os.path.join(folder, name)
    else:
        base = re.sub(r"\.(gz|old)$", "", folder)
    parts = [p for p in (base + ".gz", base + ".old", base) if os.path.isfile(p)]
    if not parts:
        raise FileNotFoundError(f"No {base}(.gz/.old) found")
    return parts


def _openText(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _hasData(path: str) -> bool:
    """False for a file that only holds the header (a run that has just started)."""
    with _openText(path) as f:
        return any(line.strip() and not re.match(r"[#A-Za-z]", line) for line in f)

def readDataHistoryPointsMultiple(
    folders: np.ndarray,
) -> Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Reads data from the HistoryPoints files in the specified folders.

    Args:
        folders: List of paths to the folders containing the HistoryPoints files.
    Returns:
        A dictionary where each key is a folder path and the value is a tuple containing:
            - points: An array of points.
            - time: An array of time values.
            - fields: An array of fields (e.g. u, v, p) at every point (axis 0) and timestep (axis 1).
    """
    data = {}
    for folder in folders:
        points, time, fields = readDataHistoryPoints(folder)
        data[folder] = (points, time, fields)
    return data


def readDataHistoryPoints(folder: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Reads data from the HistoryPoints file in the specified folder.

    Args:
        folder: Path to the folder containing the HistoryPoints file.
    Returns:
        A tuple containing:
            - points: An array of points.
            - time: An array of time values.
            - fiedls: An array of fields (e.g. u, v, p) at every point (axis 0) and timestep (axis 1).
    """
    # file_nameHistoryPoints = "HistoryPoints.his"
    # old_file_nameHistoryPoints = f"{file_nameHistoryPoints}.old"

    # oldPath = os.path.join(folder, old_file_nameHistoryPoints)
    # newPath = os.path.join(folder, file_nameHistoryPoints)

    parts = [p for p in historyParts(folder, "HistoryPoints.his") if _hasData(p)]
    read = [(p, *_readDataHistoryPoints(p)) for p in parts]
    # the newest file defines the points; older files with other points are ignored (as with the .old scheme)
    path_new, points_new, time_new, fields_new = read[-1]
    times, fieldss = [], []
    for path, points, time, fields in read:
        if points.shape != points_new.shape:
            print(
                colors.WARNING
                + f"Warning: {path} has a different number of points than {path_new}. Ignoring {path}."
                + colors.ENDC
            )
            continue
        times.append(time)
        fieldss.append(fields)
    if len(times) > 1:
        print(colors.OKGREEN + f"Info: merged {len(times)} files: {', '.join(parts)}" + colors.ENDC)

    return points_new, np.concatenate(times), np.concatenate(fieldss, axis=1)


def _readDataHistoryPoints(path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Reads data from the HistoryPoints file in the specified folder.

    Args:
        folder: Path to the folder containing the HistoryPoints file.
    Returns:
        A tuple containing:
            - points: An array of points.
            - time: An array of time values.
            - fields: An array of fields (e.g. u, v, p) at every point (axis 0) and timestep (axis 1).
    """
    # Extract data after the header
    # data: List[List[float]] = []
    # points_locations: List[List[float]] = []
    # num_points = -1
    # with open(path, "r") as f:
    #     for line in f:
    #         if not line.startswith("#"):
    #             try:
    #                 data.append(np.fromstring(line, sep=" "))
    #                 # data.append([float(value) for value in line.split()])
    #             except ValueError:
    #                 continue
    #         else:
    #             num_points += 1
    #             if num_points > 0:
    #                 points_locations.append(np.fromstring(line[1:], sep=" "))
    #                 # points_locations.append(
    #                 #     [float(value) for value in line[1:].split()][1:]
    #                 # )
    # new_data = np.array([data[i::num_points] for i in range(num_points)])

    with _openText(path) as f:
        points_locations = []
        num_points = 0

        for line in f:
            if not line.startswith("#"):
                break  # header only (a .gz holds it once, at the top)
            if line.startswith("#"):
                if num_points > 0:
                    points_locations.append(np.fromstring(line[1:], sep=" ")[1:])
                num_points += 1

    num_points -= 1

    data = np.loadtxt(path, comments="#", ndmin=2)

    new_data = data.reshape(-1, num_points, data.shape[1]).transpose(1, 0, 2)

    print(colors.OKBLUE + f"Read {len(data)} data points from {path}" + colors.ENDC)

    time = new_data[0, :, 0]  # we take the first column (time) from the first point
    variables = new_data[:, :, 1:]  # we take all columns except the first one (time)

    return np.array(points_locations), time, variables


def readErrorHistoryMultiple(
    folders: np.ndarray,
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """
    Reads energy history data from the ErrorHistory files in the specified folders.

    Args:
        folders: List of paths to the folders containing the ErrorHistory files.
    Returns:
        A dictionary where each key is a folder path and the value is a tuple containing:
            - time: An array of time values.
            - Error: An array of energy values.
    """
    data = {}
    for folder in folders:
        time, errors = readErrorHistory(folder)
        data[folder] = (time, errors)
    return data


def readErrorHistory(folder: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Reads energy history data from a file.

    Args:
        filePath: Path to the energy history file.
    Returns:
        A tuple containing:
            - time: An array of time values.
            - Error: An array of energy values.
    """
    # file_nameErrorHistory = "EnergyError.err"
    # old_file_nameErrorHistory = f"{file_nameErrorHistory}.old"
    #
    # oldPath = os.path.join(folder, old_file_nameErrorHistory)
    # newPath = os.path.join(folder, file_nameErrorHistory)

    parts = [p for p in historyParts(folder, "EnergyError.err") if _hasData(p)]
    read = [_readErrorHistory(p) for p in parts]
    time = np.concatenate([r[0] for r in read])
    errors = np.concatenate([r[1] for r in read], axis=0)

    return time, errors


def _readErrorHistory(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Reads energy history data from a file.

    Args:
        filePath: Path to the energy history file.
    Returns:
        A tuple containing:
            - time: An array of time values.
            - Error: An array of energy values.
    """
    data = np.loadtxt(path, comments=("#", "Time"), ndmin=2)  # header "Time u_L2 ..." (once per file)
    time = data[:, 0]
    errors = data[
        :, 1:
    ]  # u_L2, u_inf, u_H1, v_L2, v_inf, v_H1, p_L2, p_inf, p_H1 (9 columns)
    return time, errors


def readFourierEnergyMultiple(
    folders: np.ndarray,
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """
    Reads Fourier energy data from the EnergyFourierModes files in the specified folders.

    Args:
        folders: List of paths to the folders containing the EnergyFourierModes files.
    Returns:
        A dictionary where each key is a folder path and the value is a tuple containing:
            - time: An array of time values.
            - Energy: An 2d array of Fourier energy values, indexed by mode.
    """
    data = {}
    for folder in folders:
        time, energy = readFourierEnergy(folder)
        data[folder] = (time, energy)
    return data


def readFourierEnergy(folder: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Reads Fourier energy data from a file.

    Args:
        filePath: Path to the Fourier energy file.
    Returns:
        A tuple containing:
            - time: An array of time values.
            - Energy: An 2d array of Fourier energy values, indexed by mode.
    """

    # file_nameFourierEnergy = "EnergyFourierModes.mdl"
    # old_file_nameFourierEnergy = f"{file_nameFourierEnergy}.old"
    #
    # oldPath = os.path.join(folder, old_file_nameFourierEnergy)
    # newPath = os.path.join(folder, file_nameFourierEnergy)

    parts = [p for p in historyParts(folder, "EnergyFourierModes.mdl") if _hasData(p)]
    read = [_readFourierEnergy(p) for p in parts]
    time = np.concatenate([r[0] for r in read])
    energy = np.concatenate([r[1] for r in read], axis=1)

    return time, energy


def _readFourierEnergy(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Reads Fourier energy data from a file.

    Args:
        filePath: Path to the Fourier energy file.
    Returns:
        A tuple containing:
            - time: An array of time values.
            - Energy: An 2d array of Fourier energy values, indexed by mode.
    """
    data = np.loadtxt(path, comments="#", ndmin=2)
    time = data[:, 0]
    modes = data[:, 1]
    energy = data[:, 2]

    max_mode = int(np.max(modes)) + 1

    time_unique = np.unique(time)
    energy_array = np.reshape(energy, (len(time_unique), max_mode))
    return time_unique, energy_array.T


def editFile(
    pathFile: str, replacement_line: np.ndarray, line_startswith: np.ndarray
) -> None:
    """
    This function edits a file by replacing the lines that start with specified strings with new lines, specified in replacement_line.
    """
    # Read the file and replace the specified line
    with open(pathFile, "r") as file:
        lines = file.readlines()

    with open(pathFile, "w") as file:
        for line in lines:
            linechanged = False
            # loop in the i,j in (replacement_line, line_startswith)
            for rep_line, line_start in zip(replacement_line, line_startswith):
                if line.strip().startswith(line_start):
                    # take the last part of the line (starting with # ) and add it to rep_line
                    rep_line = rep_line + " #" + line.split("#")[-1]
                    file.write(rep_line)
                    linechanged = True
                    break
            if not linechanged:
                file.write(line)


def readFieldsBySection(dataFile: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Reads data from a file and returns three arrays:
        1. the array corresponding to the xvalues
        2. the array corresponding to the yvalues
        3. numpy array containing three indices, first corresponding to the x location, second to the y location and third to the field, for example (u, v, p, uu, uv, vv).
    """
    data = []
    xvals = np.array([])
    yvals = np.array([])

    data_dF = np.loadtxt(dataFile, skiprows=3)
    data_tmp = []
    if len(xvals) == 0:
        xvals = np.unique(data_dF[:, 0])
    if len(yvals) == 0:
        yvals = np.unique(data_dF[:, 1])
        yvals.sort()
    for x in xvals:
        tmp = data_dF[data_dF[:, 0] == x]
        data_tmp.append(tmp)
    data = np.array(data_tmp)

    return xvals, yvals, data


def extractValueXML(filename: str, param: str) -> str:
    """
    Extracts the numerical value of a parameter from an XML file, specifically
    for parameters like <P> or <p> where the parameter is defined as:
    <P> ParameterName = value </P> or <p> ParameterName = value </p>.
    
    Args:
        filename: Path to the XML file.
        param: The parameter whose value needs to be extracted.
        
    Returns:
        The value of the parameter as a float.
    """
    
    # Case-insensitive regex to find <P> or <p> tags with the parameter name and value
    p_re = re.compile(r"<\s*[Pp]\s*>\s*([A-Za-z0-9_]+)\s*=\s*([^<]+?)\s*<\s*/\s*[Pp]\s*>", re.IGNORECASE)
    
    with open(filename, "r") as f:
        for line in f:
            m = p_re.search(line)
            if m:
                param_name, value = m.groups()
                if param_name == param:
                    # Return the value as float (it should be convertible to float)
                    return value.strip()
                    try:
                        return float(value.strip())
                    except ValueError:
                        raise ValueError(f"Invalid value for '{param}': {value.strip()}")
    
    raise KeyError(f"Parameter '{param}' not found in {filename}")

