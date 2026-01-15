import re
import matplotlib.pyplot as plt
from pp.inputargs import parseArgs
import os
from pp.colors import colors as termcolors  # renamed to avoid confusion


def main():
    """
    Plot CFL and CPU time values from output.txt for each folder.
    """

    args = parseArgs()

    if len(args.folders) == 0:
        print(
            termcolors.WARNING
            + "No folders provided. Please specify at least one folder containing the output.txt file."
            + termcolors.ENDC
        )
        return

    # Prepare plot
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()

    # Standard matplotlib color cycle (tab10)
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    # Loop over folders
    for i, f in enumerate(args.folders):
        output_file = os.path.join(f, "output.txt")
        color = color_cycle[i % len(color_cycle)]

        cfl_values = []
        cpu_times = []

        # Regex patterns
        cfl_pattern = re.compile(r"CFL:\s([0-9.eE+-]+)")
        cpu_pattern = re.compile(r"CPU Time:\s([0-9.eE+-]+)s")

        try:
            with open(output_file, "r") as file:
                for line in file:
                    # Extract CFL
                    cfl_match = cfl_pattern.search(line)
                    if cfl_match:
                        cfl_values.append(float(cfl_match.group(1)))

                    # Extract CPU Time
                    cpu_match = cpu_pattern.search(line)
                    if cpu_match:
                        cpu_times.append(float(cpu_match.group(1)))

        except FileNotFoundError:
            print(
                termcolors.WARNING
                + f"File not found: {output_file}. Please ensure the file exists."
                + termcolors.ENDC
            )
            continue
        # Check matching lengths (they might differ slightly) 
        print("Lengths before trimming:") 
        print(len(cfl_values), len(cpu_times))

        # Matching lengths
        n = min(len(cfl_values), len(cpu_times))
        cfl_values = cfl_values[:n]
        cpu_times = cpu_times[:n]


        # Plot CFL
        ax1.plot(
            cfl_values,
            label=f"CFL ({os.path.basename(f)})",
            color=color,
            linestyle="--",
            linewidth=1.3,
            alpha=0.5,
        )

        # Plot CPU Time
        ax2.plot(
            cpu_times,
            label=f"CPU ({os.path.basename(f)})",
            color=color,
            linewidth=1.2,
        )

    # Labels and grid
    ax1.set_xlabel("Step Index")
    ax1.set_ylabel("CFL Value")
    ax2.set_ylabel("CPU Time (s)")
    ax1.grid(True, alpha=0.4)

    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    plt.legend(lines1 + lines2, labels1 + labels2, loc="upper right")

    plt.title("CFL and CPU Time over Steps")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

