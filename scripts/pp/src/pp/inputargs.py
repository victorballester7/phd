import argparse
import numpy as np


def parseArgs() -> argparse.Namespace:
    """
    Parse command line arguments for multiple python scripts in a folder.
    Returns:
        argparse.Namespace: Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Parse arguments for multiple python scripts in a folder."
    )

    parser.add_argument(
        "folders", 
        nargs="+", 
        help="List of folders to compare."
    )

    parser.add_argument(
        "--use_DELTA_N",
        action="store_true",
        help="Use DELTA_N with reference point the first folder.",
    )

    parser.add_argument(
        "--time_min",
        default=0.0,
        type=float,
        help="Minimum time value to include in the comparison (default: None).",
    )
    parser.add_argument(
        "--time_max",
        default=np.inf,
        type=float,
        help="Maximum time value to include in the comparison (default: None).",
    )

    parser.add_argument(
        "--points",
        metavar="int",
        nargs="+",
        type=int,
        default=[0],
        help="Point number to compare (default: [0]).",
    )

    parser.add_argument(
        "--log",
        action="store_true",
        help="Use logarithmic plotting (default: False).",
    )

    parser.add_argument(
        "--growthrate",
        action="store_true",
        help="Compute the growth rate and frequency from the data (default: False).",
    )

    parser.add_argument(
        "--testgaussian",
        action="store_true",
        help="Test if the data is Gaussian (default: False).",
    )
    
    parser.add_argument(
        "--meanmode",
        action="store_true",
        help="Include the mean mode (k=0) in the comparison (default: False).",
    )

    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress showing the plots (default: False).",
    )

    # fft and dynamicalSystem are incompatible together
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--fft",
        action="store_true",
        help="Use FFT to get frequency data (default: False).",
    )
    group.add_argument(
        "--psd",
        action="store_true",
        help="Use PSD to get frequency data (default: False).",
    )

    group.add_argument(
        "--returnMap",
        action="store_true",
        help="Plot the u vs v phase space for dynamical systems and plots a return map (default: False).",
    )
    group.add_argument(
        "--fit",
        action="store_true",
        help="Fit the data to curve of the form a * exp(b * t) * cos(c * t + d) (default: False).",
    )

    return parser.parse_args()
