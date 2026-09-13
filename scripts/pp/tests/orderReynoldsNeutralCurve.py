#!/usr/bin/env python3

from collections import OrderedDict


def reorder_block(lines, mode):
    """
    mode = 1 -> high forward + low backward
    mode = 2 -> low backward + high forward
    """

    groups = OrderedDict()

    for line in lines:
        if not line.strip():
            continue
        if line.lstrip().startswith("#"):
            continue

        parts = line.split()
        Re = float(parts[0])

        groups.setdefault(Re, []).append(line)

    high = []
    low = []

    reynolds = list(groups.keys())

    for Re in reynolds:
        entries = groups[Re]

        if len(entries) == 1:
            # Only one point exists
            high.append(entries[0])
        else:
            # Sort by frequency (2nd column)
            entries = sorted(entries, key=lambda s: float(s.split()[1]))
            low.append(entries[0])
            high.append(entries[-1])

    if mode == 1:
        return high + low[::-1]
    else:
        return low[::-1] + high


def split_blocks(lines):
    """Split file into blocks separated by blank lines."""
    blocks = []
    current = []

    for line in lines:
        if line.strip() == "":
            if current:
                blocks.append(current)
                current = []
        else:
            current.append(line)

    if current:
        blocks.append(current)

    return blocks


def reorder_continuous(lines):
    groups = OrderedDict()

    header = []

    for line in lines:
        if not line.strip():
            continue

        if line.lstrip().startswith("#"):
            header.append(line)
            continue

        parts = line.split()
        Re = float(parts[0])

        groups.setdefault(Re, []).append(line)

    low = []
    high = []

    # Reynolds are already encountered in increasing order
    for Re, entries in groups.items():
        if len(entries) == 1:
            # Only one frequency exists
            high.append(entries[0])
        else:
            # Sort by frequency
            entries = sorted(entries, key=lambda s: float(s.split()[1]))
            low.append(entries[0])
            high.append(entries[1])

    return header + low + ["\n"] + high


def main(input_file, output_file, upperLowerBranch):

    with open(input_file) as f:
        lines = f.readlines()

    # --------------------------------------------------
    # Case 1: convert a continuous loop into
    # lower branch + upper branch
    # --------------------------------------------------
    if upperLowerBranch:
        reordered = reorder_continuous(lines)

        with open(output_file, "w") as f:
            f.writelines(reordered)

        return

    # --------------------------------------------------
    # Case 2: convert two separate branches into
    # clockwise ordering
    # --------------------------------------------------

    # Keep header lines
    header = []
    idx = 0
    while idx < len(lines):
        if lines[idx].lstrip().startswith("#"):
            header.append(lines[idx])
            idx += 1
        else:
            break

    blocks = split_blocks(lines[idx:])

    if len(blocks) != 2:
        raise RuntimeError(f"Expected exactly 2 data blocks, found {len(blocks)}.")

    block1 = reorder_block(blocks[0], mode=1)
    block2 = reorder_block(blocks[1], mode=2)

    with open(output_file, "w") as f:
        f.writelines(header)

        for line in block1:
            f.write(line)

        f.write("\n")

        for line in block2:
            f.write(line)


if __name__ == "__main__":
    # cases = ["flat"]
    # cases = ["d1.5_w10", "d1.5_w15", "d1.5_w20", "flat"]
    cases = ["bfs"]
    upperLowerBranch = True
    leftRightBranch = not upperLowerBranch

    for c in cases:
        input_file = f"../../data/neutralCurves/{c}/neutral_curve_temporal_gaster.dat"
        output_file = input_file

        if upperLowerBranch:
            main(input_file, output_file, True)
        else:
            main(input_file, output_file, False)
