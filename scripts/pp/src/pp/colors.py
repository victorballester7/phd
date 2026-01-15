import numpy as np
import matplotlib.colors as mcolors

class colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def mix_with_black(color, alpha=0.5):
    """
    Mix a color with black using alpha blending.
    
    Parameters
    ----------
    color : str or tuple
        Color to mix with black
    alpha : float, optional
        Blending factor (0 = black, 1 = original color)
        
    Returns
    -------
    tuple
        RGB tuple of mixed color
    """
    rgb = np.array(mcolors.to_rgb(color))
    return (1 - alpha) * rgb
