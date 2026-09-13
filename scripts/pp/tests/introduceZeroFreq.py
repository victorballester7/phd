#!/usr/bin/env python3

import os
import re
import sys

from pp.fileManagement import extract_depth_width


def main(filename):

    _, width = extract_depth_width(filename)

    # Read file
    rows = []
    header_line = ""
    with open(filename) as f:
        for line in f:
            if not line.strip():
                continue
            if line.lstrip().startswith("#"):
                header_line = line
                continue
            Re, freq, freq2, x = map(float, line.split())
            rows.append([Re, freq, x])

    new_rows = []


    # All x locations inside the gap
    x_values = sorted(set(r[2] for r in rows if 0 <= r[2] <= width))

    for x_target in x_values:
        matches = [r for r in rows if abs(r[2] - x_target) < 1e-8]

        Re = matches[0][0]

        if len(matches) == 1:
            new_rows.append(matches[0])

            if abs(matches[0][1]) > 1e-12:
                new_rows.append([Re, 0.0, x_target])

        else:
            keep = max(matches, key=lambda r: r[1])
            new_rows.append(keep)
            new_rows.append([Re, 0.0, x_target])

    # Add rows for x outside the gap unchanged
    for row in rows:
        if row[2] < 0 or row[2] > width:
            new_rows.append(row)

    # Sort by x then frequency
    new_rows.sort(key=lambda r: (r[2], r[1]))

    with open(filename, "w") as f:
        f.write(header_line)
        for Re, freq, x in new_rows:
            f.write(f"{Re:.6f} {freq:.6f} {x:.6f}\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: uv run {sys.argv[0]} datafile")
        sys.exit(1)

    main(sys.argv[1])
