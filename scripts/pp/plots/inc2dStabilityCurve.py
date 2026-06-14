import os
import matplotlib.pyplot as plt
from pp.colors import mix_with_black, colors
from pp.analyzeData import shared_depth_curve, build_bifurcation_pairs
import pandas as pd
from pp.smoothing import smooth_curve
import matplotlib as mpl
from pp.fontSizeLaTex import compute_mpl_fontsize, compute_figure_size


# plt.style.use('plots/style/tsfp.mplstyle')
plt.style.use('plots/style/customvictor.mplstyle')



MARKERS = ["o", "^", "d"]
MARKERSIZE = 10
MARKER_COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"][:4]

# labels as dictionary
LABELS = {
        "equilibrium": "Equilibrium",
        "limitcycle": "Limit cycle",
        "chaoticattractor": "Chaotic attractor",
}

BIFURCATION_LABELS = [
    "Hopf bifurcation",
    "Chaos onset",
]

LINE_STYLES = [":", "-.", "--"]


def plot_stability_diagram(df, save_path):
    """
    Create and save the stability bifurcation diagram.

    Parameters
    ----------
    data_sets : list of np.ndarray
        List of stability region datasets
    save_path : str
        Path to save the output PDF
    """
    # LaTeX document settings: 10pt font, figure width 6cm
    latex_width_cm = 7.85 # linewidth of tsfp template
    latex_font_pt = 9 # 9pt small font size in tsfp template (figures captions)
    
    fig_width_in, fig_height_in = compute_figure_size(latex_width_cm, aspect_ratio=0.75)

    # we make the figure as if there were only one column
    fig_width_in *= 2
    fig_height_in *= 2
    print(f"Figure size (inches): {fig_width_in:.2f} x {fig_height_in:.2f}")
    compute_mpl_fontsize(fig_width_in=fig_width_in, latex_width_cm=latex_width_cm, latex_font_pt=latex_font_pt)
    fig, ax = plt.subplots(figsize=(fig_width_in, fig_height_in), constrained_layout=True)

    
    inc2dData = df.loc[df["run"] == "inc2dRe1000"].to_numpy()

    for m, c, lab in zip(MARKERS, MARKER_COLORS, LABELS.keys()):
        subset = inc2dData[inc2dData[:, 4] == lab]
        edge_color = mix_with_black(c, alpha=0.3)
        ax.plot(
                subset[:, 0],
                subset[:, 1],
                marker=m,
                markersize=MARKERSIZE,
                markerfacecolor=c,
                markeredgecolor=edge_color,
                linestyle="none",
                label=LABELS[lab],
                zorder=3,
        )

    # Plot bifurcation curves
    # bifurcation_pairs = build_bifurcation_pairs(inc2dData)

    # for pair, line_style, label in zip(
    #     bifurcation_pairs, LINE_STYLES, BIFURCATION_LABELS
    # ):
    #     w_curve, d_curve = shared_depth_curve(pair[0], pair[1])
    #     if len(w_curve) > 4:
    #         w_curve, d_curve = smooth_curve(w_curve, d_curve)
    #     # Smooth the curve using spline interpolation
    #     # if len(w_curve) >= 4:  # Require at least 4 points for spline
    #     #     spline = UnivariateSpline(d_curve, w_curve, s=1.0)
    #     #     d_smooth = np.linspace(np.min(d_curve), np.max(d_curve), 300)
    #     #     w_smooth = spline(d_smooth)
    #     #     w_curve, d_curve = w_smooth, d_smooth

    #     ax.plot(
    #         w_curve,
    #         d_curve,
    #         linestyle=line_style,
    #         label=label,
    #         zorder=2,
    #     )




    # Axis labels and formatting
    ax.set_xlabel(r"$w/\delta^*$")
    ax.set_ylabel(r"$d/\delta^*$", rotation=0, labelpad=20)

    # Grid and legend
    ax.legend()

    # Tight layout
    plt.savefig(save_path, format="pdf")
    print(colors.OKGREEN + f"✓ Plot saved to: {save_path}" + colors.ENDC)

def main():
    """Main execution function."""
    # Set up paths
    script_path = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_path, "../../../data/stability/stability.dat")
    output_path = os.path.join(script_path, "../../../images/inc2dStabilityCurve.pdf")

    # Load data
    data_frames = pd.read_csv(data_dir, sep=" ", comment="#", names=["w", "d", "sigma", "omega", "type", "run"])

    # Create plot
    plot_stability_diagram(data_frames, output_path)


if __name__ == "__main__":
    main()
