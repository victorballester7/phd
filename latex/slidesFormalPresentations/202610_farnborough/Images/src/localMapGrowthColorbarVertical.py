"""Vertical version of localMapGrowthColorbar.pdf for the slides.

Colour map, levels and ticks are imported from the paper's plotting script, so
the bar stays in sync with the local stability maps. Run from scripts/pp:

    uv run python ~/Desktop/PhD/latex/slidesFormalPresentations/202610_farnborough/Images/src/localMapGrowthColorbarVertical.py
"""

import os

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import FixedFormatter

import plots.localNfactorMaps as lm

W_CM, H_CM = 1.35, 3.6
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../localMapGrowthColorbarVertical.pdf")

fig, ax = plt.subplots(figsize=(W_CM / 2.54, H_CM / 2.54), layout="constrained")
fig.get_layout_engine().set(w_pad=0.015, h_pad=0.02, wspace=0, hspace=0)
ax.grid(False)
cb = fig.colorbar(
    ScalarMappable(norm=lm.GROWTH_NORM, cmap=lm.GROWTH_CMAP),
    cax=ax,
    orientation="vertical",
    spacing="uniform",
    ticks=lm.GROWTH_TICKS,
)
cb.set_label(r"$-\alpha_\mathrm{i}$", rotation=0, va="center", ha="left", labelpad=0.5)
cb.ax.yaxis.set_major_formatter(FixedFormatter(lm.GROWTH_TICKLABELS))
cb.ax.tick_params(which="both", direction="out", length=2.0, width=0.5, pad=1.5)
cb.outline.set_linewidth(0.5)
cb.dividers.set_visible(False)
fig.savefig(OUT, format="pdf")
print("Saved", os.path.normpath(OUT))
