from pathlib import Path

TIME_SHIFT = 6720

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def process_file(fname):
    fname = Path(fname)
    outname = fname.with_suffix(fname.suffix + ".new")

    with open(fname, "r") as fin, open(outname, "w") as fout:
        for line in fin:
            stripped = line.strip()

            # Preserve empty lines and comments
            if not stripped or stripped.startswith("#"):
                fout.write(line)
                continue

            parts = line.split()

            # Case 1: EnergyError.err or HistoryPoints.his
            # → time is first column and is float
            if is_float(parts[0]):
                t = float(parts[0]) + TIME_SHIFT
                parts[0] = f"{t:.16e}" if "e" in parts[0].lower() else f"{t:.6f}"
                fout.write(" ".join(parts) + "\n")
                continue

            # Case 2: EnergyFourierModes.mdl
            # → format: time, mode, energy
            if len(parts) >= 3 and is_float(parts[0]):
                t = float(parts[0]) + TIME_SHIFT
                parts[0] = f"{t:g}"
                fout.write(" ".join(parts) + "\n")
                continue

            # Fallback: copy line
            fout.write(line)

    print(f"Wrote {outname}")

if __name__ == "__main__":
    import sys
    for f in sys.argv[1:]:
        process_file(f)

