import argparse
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "../../../data/psd"
IMAGE_DIR = SCRIPT_DIR / "../../../images"
FILE_RE = re.compile(r"d2_w(?P<width>[0-9]+(?:\.[0-9]+)?)_Re800_x[0-9]+(?:\.[0-9]+)?_y[0-9]+(?:\.[0-9]+)?_z[0-9]+(?:\.[0-9]+)?\.dat$")
DEFAULT_VARIABLE = "v"


def _style_path():
    style = SCRIPT_DIR / "style/customvictor.mplstyle"
    return str(style) if style.exists() else "default"


def _load_psd_file(path, variable):
    header = path.read_text().splitlines()[0].lstrip("#").split()
    data = np.loadtxt(path)

    matches = [
        idx
        for idx, name in enumerate(header)
        if name == f"PSD_{variable}" or name.startswith(f"PSD_{variable}_")
    ]
    if not matches:
        available = ", ".join(name for name in header if name.startswith("PSD_"))
        raise ValueError(
            f"Could not find PSD column for variable '{variable}' in {path.name}. "
            f"Available PSD columns: {available}"
        )

    omega = data[:, 0]
    psd = data[:, matches[0]]
    return omega, psd


def _width_edges(widths):
    widths = np.asarray(widths, dtype=float)
    if widths.size == 1:
        dx = max(abs(widths[0]) * 0.05, 1.0)
        return np.array([widths[0] - dx / 2, widths[0] + dx / 2])

    edges = np.empty(widths.size + 1)
    edges[1:-1] = 0.5 * (widths[:-1] + widths[1:])
    edges[0] = widths[0] - 0.5 * (widths[1] - widths[0])
    edges[-1] = widths[-1] + 0.5 * (widths[-1] - widths[-2])
    return edges


def _omega_edges(omega):
    omega = np.asarray(omega, dtype=float)
    if omega.size == 1:
        return np.array([omega[0] - 0.5, omega[0] + 0.5])

    edges = np.empty(omega.size + 1)
    edges[1:-1] = 0.5 * (omega[:-1] + omega[1:])
    edges[0] = omega[0] - 0.5 * (omega[1] - omega[0])
    edges[-1] = omega[-1] + 0.5 * (omega[-1] - omega[-2])
    return edges


def load_psd_by_width(data_dir, variable):
    records = []
    for path in sorted(data_dir.glob("d2_w*_Re800_*.dat")):
        match = FILE_RE.match(path.name)
        if match is None:
            continue

        width = float(match.group("width"))
        omega, psd = _load_psd_file(path, variable)
        records.append((width, omega, psd, path))

    if not records:
        raise FileNotFoundError(f"No d2_w*_Re800.dat files found in {data_dir}")

    records.sort(key=lambda record: record[0])
    omega_ref = records[0][1]
    for _, omega, _, path in records[1:]:
        if omega.shape != omega_ref.shape or not np.allclose(omega, omega_ref):
            raise ValueError(
                f"Omega grid in {path.name} does not match {records[0][3].name}"
            )

    widths = np.array([record[0] for record in records])
    amplitudes = np.column_stack([record[2] for record in records])
    max_amplitude = np.nanmax(amplitudes)
    if not np.isfinite(max_amplitude) or max_amplitude <= 0:
        raise ValueError("The selected PSD amplitudes have no positive finite values")

    return widths, omega_ref, amplitudes / max_amplitude


def make_plot(variable=DEFAULT_VARIABLE, data_dir=DATA_DIR, output=None):
    data_dir = Path(data_dir).resolve()
    widths, omega, psd = load_psd_by_width(data_dir, variable)
    # psd shape: (n_omega, n_widths)

    # ── Per-slice normalisation ───────────────────────────────────────────
    col_max = np.nanmax(psd, axis=0, keepdims=True)   # (1, n_widths)
    col_max[col_max == 0] = 1.0                        # guard against silent slices
    psd_norm = psd / col_max                           # each column in [0, 1]

    # ── Log scale: reveal sub-harmonics ──────────────────────────────────
    eps = 1e-6                                         # floor to avoid log(0)
    psd_db = np.log(np.clip(psd_norm, eps, 1)) # in dB, range [-60, 0]
    # psd_db = 10 * np.log10(np.clip(psd_norm, eps, 1)) # in dB, range [-60, 0]

    plt.style.use(_style_path())
    fig, ax = plt.subplots(constrained_layout=True)

    # ── Plot each w slice as a narrow independent strip ──────────────────────────
    STRIP_WIDTH = 2.0   # in data units (w/δ*); tweak to taste

    for i, w in enumerate(widths):
        w_edges = np.array([w - STRIP_WIDTH / 2, w + STRIP_WIDTH / 2])
        mesh = ax.pcolormesh(
            w_edges,
            _omega_edges(omega),
            psd_db[:, i : i + 1],   # (n_omega, 1) — single column
            cmap="inferno",
            # vmin=-40,
            # vmax=0,
            shading="flat",
            rasterized=True,
        )
    # mesh = ax.pcolormesh(
    #     _width_edges(widths),
    #     _omega_edges(omega),
    #     psd_db,
    #     cmap="inferno",       # perceptually uniform, prints well in B&W
    #     vmin=-40,             # cut off noise floor at -40 dB
    #     vmax=0,
    #     shading="flat",
    #     rasterized=True,      # essential for PDF — avoids huge vector files
    # )

    ax.set_facecolor("black")

    ax.set_ylim([0, 1])
    ax.set_xlabel(r"$w/\delta^*$")
    ax.set_ylabel(r"$\omega$", rotation=0, labelpad=10)

    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label(
        rf"$\log\!\left(S_{{{variable}}}/\max_\omega S_{{{variable}}}\right)$"
    )


    if output is None:
        output = IMAGE_DIR / f"psd_vs_width_{variable}.pdf"
    else:
        output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, format=output.suffix.lstrip(".") or "pdf")
    print(f"Plot saved to: {output}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot omega vs width, colored by normalized PSD amplitude."
    )
    parser.add_argument(
        "-v",
        "--variable",
        default=DEFAULT_VARIABLE,
        help="PSD variable to use: u, v, p, ... (default: v)",
    )
    parser.add_argument(
        "--data-dir",
        default=DATA_DIR,
        type=Path,
        help="Directory containing d2_w*_Re800.dat files",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output figure path (default: ../../../images/psd_vs_width_<variable>.pdf)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    make_plot(variable=args.variable, data_dir=args.data_dir, output=args.output)
