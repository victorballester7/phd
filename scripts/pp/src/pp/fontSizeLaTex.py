import matplotlib as mpl


def compute_figure_size(latex_width_cm, aspect_ratio=0.75):
    """
    Compute figure size in inches matching the LaTeX target width.

    Using the same physical size avoids extreme scaling and keeps
    font rendering crisp.

    Parameters
    ----------
    latex_width_cm : float
        Width in LaTeX document (cm)
    aspect_ratio : float
        height / width ratio (default 0.75 = 4:3)

    Returns
    -------
    tuple : (width_in, height_in)
    """
    cm_to_inch = 1 / 2.54
    width_in = latex_width_cm * cm_to_inch
    height_in = width_in * aspect_ratio
    return width_in, height_in


def compute_mpl_fontsize(latex_width_cm, latex_font_pt):
    """
    Compute matplotlib font size so that, after scaling in LaTeX,
    it matches the LaTeX document font size.

    When LaTeX scales the figure, text scales proportionally. If the figure
    is shrunk by factor S, the text also shrinks by S. To get final text at
    latex_font_pt, we need: mpl_fontsize * S = latex_font_pt
    => mpl_fontsize = latex_font_pt / S

    Parameters
    ----------
    fig_width_in : float
        Width of figure in inches (e.g. 8 from figsize=(8,6))
    latex_width_cm : float
        Width you will use in LaTeX (e.g. 6 cm)
    latex_font_pt : float
        Document font size (e.g. 10, 11, 12)

    Returns
    -------
    float : matplotlib font size to use
    """

    cm_to_inch = 1 / 2.54
    latex_width_in = latex_width_cm * cm_to_inch

    fig_width_in = mpl.rcParams["figure.figsize"][0]  # width in inches

    # S = how much LaTeX will scale the figure (< 1 means shrink)
    ratio_widthFontsize_mpl = fig_width_in / mpl.rcParams["font.size"]
    ratio_widthFontsize_latex = latex_width_in / latex_font_pt

    print(f"ratio_widthFontsize_mpl: {ratio_widthFontsize_mpl:.3f}")
    print(f"ratio_widthFontsize_latex: {ratio_widthFontsize_latex:.3f}")

    scale_factor = ratio_widthFontsize_mpl / ratio_widthFontsize_latex
    # scale_factor = 1/scale_factor

    # To get latex_font_pt after scaling: mpl_fontsize * scale_factor = latex_font_pt
    # scale_factor = fig_width_in / latex_width_in

    print(
        f"Scale factor: {scale_factor:.3f}"
    )

    mpl.rcParams.update(
        {
            "font.size": mpl.rcParams["font.size"] * scale_factor,
            "axes.labelsize": mpl.rcParams["axes.labelsize"] * scale_factor,
            "axes.titlesize": mpl.rcParams["axes.titlesize"] * scale_factor,
            "axes.linewidth": mpl.rcParams["axes.linewidth"] * scale_factor,
            "grid.linewidth": mpl.rcParams["grid.linewidth"] * scale_factor,
            "lines.linewidth": mpl.rcParams["lines.linewidth"] * scale_factor,
            "lines.markersize": mpl.rcParams["lines.markersize"] * scale_factor,
            "lines.markeredgewidth": mpl.rcParams["lines.markeredgewidth"] * scale_factor,
            "legend.fontsize": mpl.rcParams["legend.fontsize"] * scale_factor,
            "xtick.labelsize": mpl.rcParams["xtick.labelsize"] * scale_factor,
            "ytick.labelsize": mpl.rcParams["ytick.labelsize"] * scale_factor,
            "xtick.major.width": mpl.rcParams["xtick.major.width"] * scale_factor,
            "ytick.major.width": mpl.rcParams["ytick.major.width"] * scale_factor,
        }
    )

    print("axes.grid =", mpl.rcParams["axes.grid"])
