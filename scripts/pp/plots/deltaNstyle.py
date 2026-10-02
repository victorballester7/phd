"""
Shared styling for the Delta n contour maps (Re = 800, 1000, 3000).

The three maps live in two different scripts, so the colour scale and the
treatment of the globally unstable region are defined here once. They have to
agree: the maps are read against each other, and the reader has no way of
knowing that two shades of the same colour mean the same number unless they
come from the same ramp.

Two choices are deliberate.

*Dark means large.* The earlier version used a reversed ramp, so the largest
amplification came out the palest. That is backwards, and it collided with the
second problem below.

*The excluded region is grey, not white.* Above the Hopf boundary there is no
steady base flow and no Delta n to plot. Filling that region with white made it
identical to the top of a reversed colour scale, so "no data" and "strongest
amplification" looked the same. A flat grey is not a value on the ramp and
cannot be mistaken for one.
"""

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

from pp.filterData import average_curves
from pp.hopfDataPoints import hopf2d_re_stable, hopf_re_unstable

# One ramp for every Reynolds number. The top is truncated well short of full
# saturation so that the contour lines drawn on top stay legible at the deep
# end of the scale.
BASE_CMAP = "Blues"
CMAP_RANGE = (0.03, 0.82)

COL_NUMERIC = "black"
COL_EXPERIMENT = "tab:orange"
COL_BOUNDARY = "black"
COL_EXCLUDED = "0.9"

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=256):
    return mcolors.LinearSegmentedColormap.from_list(
        f"trunc({cmap.name},{minval:.2f},{maxval:.2f})",
        cmap(np.linspace(minval, maxval, n)),
    )


def deltaN_colormap():
    return truncate_colormap(plt.get_cmap(BASE_CMAP), *CMAP_RANGE)


def add_numeric_contours(ax, xgrid, ygrid, zgrid, levels, label_positions):
    cs = ax.contour(
        xgrid,
        ygrid,
        zgrid,
        levels=levels,
        colors=COL_NUMERIC,
        linewidths=0.8,
        zorder=3,
    )
    ax.clabel(cs, inline=True, manual=label_positions, fmt="%d")
    return cs


def add_excluded_region(ax, re, label=None):
    """
    Hatch the part of the (w, d) plane with no steady base flow.

    The shading is drawn *over* the contour lines, as the white mask it
    replaces was: a level set of Delta n has no meaning where there is no
    steady base flow to compute it from.

    ``label`` is passed to the boundary line so the caller can decide whether
    the panel carries a legend.
    """
    bif = average_curves(hopf2d_re_stable[re], hopf_re_unstable[re], npts=500)[2]
    w, d = bif[:, 0], bif[:, 1]

    ax.fill_between(
        w,
        d,
        y2=ax.get_ylim()[1],
        facecolor=COL_EXCLUDED,
        edgecolor="none",
        linewidth=0.0,
        zorder=5,
    )
    ax.plot(
        w,
        d,
        color=COL_BOUNDARY,
        linestyle=":",
        linewidth=0.9,
        label=label,
        zorder=6,
    )
