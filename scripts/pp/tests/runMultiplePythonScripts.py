import os
import subprocess


def main():
    script = "scripts/plotHistoryPoints.py"
    mainFolder = "~/hosts/hpc/PhD/runs/incGapRe800/baseflow/dns/"

    # folder, point
    cases = [
        ["d2_w44", 27],
        ["d2_w47", 29],
        ["d2_w53", 31],
        ["d2_w57", 33],
        ["d2_w64", 35],
        ["d2_w70", 39],
        ["d2_w72", 39],
        ["d2_w75", 41],
        ["d2_w78", 41],
        ["d2_w80", 43],
        ["d2_w82.5", 43],
        ["d2_w85", 45],
        ["d2_w90", 47],
        ["d2_w95", 49],
        ["d2_w98", 49],
        ["d2_w102", 51],
        ["d2_w105", 53],
        ["d2_w108", 53],
        ["d2_w114", 55],
        ["d2_w120", 59],
        ["d2_w125", 61],
    ]

    for f, p in cases:
        folder = os.path.join(mainFolder, f)
        command = f"uv run {script} {folder} --point {p} --time_min 2000 --psd --quiet"
        print(f"Running command: {command}")
        subprocess.run(command, shell=True, check=True)


if __name__ == "__main__":
    main()
