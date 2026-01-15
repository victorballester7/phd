import os
import matplotlib.pyplot as plt
from pp.colors import mix_with_black, colors
from pp.analyzeData import shared_depth_curve, build_bifurcation_pairs
import pandas as pd
from pp.smoothing import smooth_curve

plt.style.use('plots/style/customvictor.mplstyle')

# Data and plot configuration
DATA_FILES = [
    "01_equilibrium.dat",
    "02_po.dat",
    "03_pocascade.dat",
    "04_chaos.dat",
]


MARKERS = ["o", "^", "d", "s"]
MARKER_COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"][:4]

LABELS = [
    "Equilibrium",
    "P.O. of fundamental period",
    "Period-doubling cascade",
    "Chaos",
]

BIFURCATION_LABELS = [
    "Hopf bifurcation",
    "Period-doubling bifurcation",
    "Chaos onset",
]

LINE_STYLES = [":", "-.", "--"]


def plot_stability_diagram(data_frames, save_path):
    """
    Create and save the stability bifurcation diagram.

    Parameters
    ----------
    data_sets : list of np.ndarray
        List of stability region datasets
    save_path : str
        Path to save the output PDF
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    inc2dData = [
        df.loc[df["run"] == "inc2d", ["w", "d"]].to_numpy() for df in data_frames
    ]
    
    # Plot stability regions with markers
    for dataset, marker, color, label in zip(inc2dData, MARKERS, MARKER_COLORS, LABELS):
        edge_color = mix_with_black(color, alpha=0.3)
        ax.plot(
            dataset[:, 0],
            dataset[:, 1],
            marker=marker,
            markerfacecolor=color,
            markeredgecolor=edge_color,
            linestyle="none",
            label=label,
            zorder=3,
        )

    # Plot bifurcation curves
    bifurcation_pairs = build_bifurcation_pairs(inc2dData)

    for pair, line_style, label in zip(
        bifurcation_pairs, LINE_STYLES, BIFURCATION_LABELS
    ):
        w_curve, d_curve = shared_depth_curve(pair[0], pair[1])
        if len(w_curve) > 4:
            w_curve, d_curve = smooth_curve(w_curve, d_curve)
        # Smooth the curve using spline interpolation
        # if len(w_curve) >= 4:  # Require at least 4 points for spline
        #     spline = UnivariateSpline(d_curve, w_curve, s=1.0)
        #     d_smooth = np.linspace(np.min(d_curve), np.max(d_curve), 300)
        #     w_smooth = spline(d_smooth)
        #     w_curve, d_curve = w_smooth, d_smooth

        ax.plot(
            w_curve,
            d_curve,
            linestyle=line_style,
            label=label,
            zorder=2,
        )




    # Axis labels and formatting
    ax.set_xlabel(r"$w/\delta^*$")
    ax.set_ylabel(r"$d/\delta^*$", rotation=0, labelpad=20)

    # Grid and legend
    ax.legend()

    # Tight layout and save
    plt.savefig(save_path)
    print(colors.OKGREEN + f"✓ Plot saved to: {save_path}" + colors.ENDC)

def main():
    """Main execution function."""
    # Set up paths
    script_path = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_path, "../../../data/stability/")
    output_path = os.path.join(script_path, "../../../images/inc2dStabilityCurve.pdf")

    # Load data
    file_paths = [os.path.join(data_dir, f) for f in DATA_FILES]
    data_frames = [
        pd.read_csv(f, sep=" ", comment="#", names=["w", "d", "sigma", "omega", "run"])
        for f in file_paths
    ]

    # Create plot
    plot_stability_diagram(data_frames, output_path)


if __name__ == "__main__":
    main()
