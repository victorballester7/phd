import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from pp.colors import colors

CM_TO_INCH = 1 / 2.54
MAX_WIDTH_JFM = 13.5  # max width an image can be (cm)


def jfm_ticks(ax, n_minor=2, top=True):
    """
    Minor ticks between the major ones and ticks mirrored on all four sides.

    top=False leaves the top edge alone, for an axes whose top carries a
    secondary x axis with its own ticks.

    Call it after the axis scales are set: a log axis keeps its own minor
    ticks instead of getting linear ones.
    """
    for axis in (ax.xaxis, ax.yaxis):
        if axis.get_scale() == "linear":
            axis.set_minor_locator(AutoMinorLocator(n_minor))
    ax.tick_params(which="both", top=top, right=True)
    ax.tick_params(which="minor", length=2)


def axis_on_top(ax, zorder=10):
    """Draw ticks and spines above every artist in the axes."""
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_zorder(zorder)
    for spine in ax.spines.values():
        spine.set_zorder(zorder)


class FigureFrame:
    """
    Figure with a fixed-size axes frame and fixed paddings around it.

    All geometry is given in cm. The figure size is built from the frame
    size plus the paddings, so the axes box has the same physical size in
    every figure regardless of labels, legends, etc.

    Parameters
    ----------
    frame_w : float
        Width of the axes box (cm)
    aspect_ratio : float
        Frame height / frame width
    pad_l, pad_r : float
        Space to the left / right of the frame (cm); pad_r defaults to pad_l
    pad_b, pad_t : float
        Space below / above the frame (cm)
    max_width : float
        Maximum allowed figure width (cm); a warning is printed if exceeded

    Attributes
    ----------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    """

    def __init__(
        self,
        frame_w=6.0,
        aspect_ratio=0.75,
        pad_l=0.0,
        pad_r=None,
        pad_b=0.0,
        pad_t=0.0,
        max_width=MAX_WIDTH_JFM,
    ):
        if pad_r is None:
            pad_r = pad_l
        frame_h = frame_w * aspect_ratio
        fig_w = frame_w + pad_l + pad_r
        fig_h = frame_h + pad_b + pad_t

        if fig_w > max_width:
            print(
                colors.WARNING
                + f"Warning: Figure width {fig_w:.2f} cm exceeds {max_width} cm, which may not fit in a single column of a journal."
                + colors.ENDC
            )
        else:
            print(
                colors.OKGREEN
                + f"Figure width {fig_w:.2f} cm is within the limit of {max_width} cm."
                + colors.ENDC
            )

        self.frame_w, self.frame_h = frame_w, frame_h
        self.pad_l, self.pad_b = pad_l, pad_b
        self.fig_w, self.fig_h = fig_w, fig_h

        self.fig = plt.figure(figsize=(fig_w * CM_TO_INCH, fig_h * CM_TO_INCH))
        self.ax = self.fig.add_axes(
            [pad_l / fig_w, pad_b / fig_h, frame_w / fig_w, frame_h / fig_h]
        )

    def jfm_ticks(self, n_minor=2, top=True):
        """Apply jfm_ticks to the frame axes."""
        jfm_ticks(self.ax, n_minor, top)

    def add_axes_cm(self, left, bottom, width, height):
        """Add an extra axes placed in cm from the lower-left figure corner."""
        return self.fig.add_axes(
            [
                left / self.fig_w,
                bottom / self.fig_h,
                width / self.fig_w,
                height / self.fig_h,
            ]
        )

    def axis_on_top(self, zorder=10):
        """Apply axis_on_top to the frame axes."""
        axis_on_top(self.ax, zorder)

    def colorbar(self, mappable, width=0.3, gap=0.2, **kwargs):
        """
        Add a colorbar to the right of the frame, inside pad_r.

        Unlike fig.colorbar(mappable, ax=ax), the frame is not shrunk:
        the colorbar gets its own axes, so pad_r must be large enough to
        hold gap + width + ticks + label.

        Parameters
        ----------
        mappable : matplotlib.cm.ScalarMappable
        width : float
            Width of the colorbar (cm)
        gap : float
            Space between the frame and the colorbar (cm)
        **kwargs
            Passed to fig.colorbar

        Returns
        -------
        matplotlib.colorbar.Colorbar
        """
        x0 = self.pad_l + self.frame_w + gap
        cax = self.add_axes_cm(x0, self.pad_b, width, self.frame_h)
        cbar = self.fig.colorbar(mappable, cax=cax, **kwargs)
        # the style's axes.grid would otherwise draw grid lines inside the bar
        cax.grid(False)
        return cbar
