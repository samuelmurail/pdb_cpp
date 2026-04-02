#!/usr/bin/env python3
# coding: utf-8

"""Generate a publication-quality benchmark figure.

Top row  – Three panels (small / medium / big protein) showing grouped bar
           charts of common operations, one bar per library, log-scale y-axis.
Bottom   – DockQ end-to-end panel, pdb_cpp vs DockQ CLI, with speedup
           annotations.

Usage
-----
PYTHONPATH=src python benchmark/plot_manuscript_figure.py
PYTHONPATH=src python benchmark/plot_manuscript_figure.py \
    --common benchmark/benchmark_common.csv \
    --dockq  benchmark/benchmark_dockq.csv \
    --output benchmark/benchmark_figure.pdf
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Colorblind-safe palette (Okabe-Ito subset)
COLORS = {
    "pdb_cpp":   "#0072B2",
    "pdb_numpy": "#56B4E9",
    "biopython": "#E69F00",
    "biotite":   "#009E73",
    "DockQ":     "#D55E00",
}

PREFERRED_LIB_ORDER = ["pdb_cpp", "pdb_numpy", "biopython", "biotite"]

# Nicer operation labels for the x-axis
OP_LABELS = {
    "read":                    "read",
    "write":                   "write",
    "select_within10_chainA":  "select\n(within 10 Å)",
    "select_complex":          "select",
    "get_aa_seq":              "sequence",
    "rmsd_ca_shift":           "RMSD",
    "dihedral_ca":             "dihedral",
    "align_seq_chainA":        "seq. align",
    "align_ca_self":           "struct.\nalign",
    "hbond":                   "H-bonds",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Create manuscript benchmark figure")
    p.add_argument(
        "--common",
        default="benchmark/benchmark_common.csv",
        help="CSV from compare_common_speed.py",
    )
    p.add_argument(
        "--dockq",
        default="benchmark/benchmark_dockq.csv",
        help="CSV from compare_dockq_speed.py",
    )
    p.add_argument(
        "--output",
        default="benchmark/benchmark_figure.pdf",
        help="Output figure path (PDF or PNG)",
    )
    return p.parse_args()


def _mean_sem(values: list[float]) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    if len(values) == 1:
        return values[0], 0.0
    arr = np.asarray(values, dtype=float)
    return float(np.mean(arr)), float(np.std(arr, ddof=1) / math.sqrt(len(arr)))


def _atom_count(path: Path) -> int:
    """Count ATOM/HETATM lines in a PDB/mmCIF file."""
    if not path.exists():
        return 0
    suffix = path.suffix.lower()
    count = 0
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if suffix == ".pdb":
                if line.startswith(("ATOM", "HETATM")):
                    count += 1
            elif suffix in (".cif", ".mmcif"):
                stripped = line.lstrip()
                if stripped.startswith(("ATOM ", "HETATM ")):
                    count += 1
    return count


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

def load_common(csv_path: Path) -> dict:
    """Return {file: {operation: {library: [mean_s, ...]}}}."""
    data: dict[str, dict[str, dict[str, list[float]]]] = {}
    with csv_path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            f = row["file"]
            op = row["operation"]
            lib = row["library"]
            data.setdefault(f, {}).setdefault(op, {}).setdefault(lib, []).append(
                float(row["mean_s"])
            )
    return data


def load_dockq(csv_path: Path) -> list[dict]:
    """Return list of row dicts from the DockQ CSV."""
    with csv_path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------------------
# Panel A — common operations
# ---------------------------------------------------------------------------

SIZE_LABELS = ["Small", "Medium", "Large"]


def _sorted_file_keys(data: dict) -> list[str]:
    """Return file keys sorted by atom count (ascending)."""
    keys = list(data.keys())
    keys.sort(key=lambda f: _atom_count(Path(f)))
    return keys


def draw_common_panel(ax: plt.Axes, data: dict, file_key: str,
                      size_label: str, show_ylabel: bool = True,
                      show_legend: bool = False) -> None:
    """Draw a grouped bar chart for one file."""
    per_op = data[file_key]
    operations = list(per_op.keys())
    libraries = sorted(
        {lib for op_d in per_op.values() for lib in op_d},
        key=lambda lib: PREFERRED_LIB_ORDER.index(lib) if lib in PREFERRED_LIB_ORDER else 99,
    )

    x = np.arange(len(operations), dtype=float)
    n_libs = max(1, len(libraries))
    width = 0.8 / n_libs

    for idx, lib in enumerate(libraries):
        means, sems = [], []
        for op in operations:
            vals = per_op.get(op, {}).get(lib, [])
            m, se = _mean_sem(vals)
            means.append(m)
            sems.append(se)
        offsets = x - 0.4 + width / 2 + idx * width
        ax.bar(
            offsets,
            means,
            width=width,
            yerr=sems,
            capsize=2,
            label=lib,
            color=COLORS.get(lib, "#888888"),
            alpha=0.9,
            edgecolor="black",
            linewidth=0.3,
        )

    pdb_id = Path(file_key).stem.upper()
    n_atoms = _atom_count(Path(file_key))
    title = f"{size_label} – {pdb_id}"
    if n_atoms:
        title += f" ({n_atoms:,} atoms)"
    ax.set_title(title, fontsize=8, fontweight="bold")
    ax.set_yscale("log")
    if show_ylabel:
        ax.set_ylabel("Execution time (s)", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [OP_LABELS.get(op, op) for op in operations],
        rotation=45,
        ha="right",
        fontsize=5.5,
    )
    ax.tick_params(axis="y", labelsize=6)
    ax.grid(True, axis="y", alpha=0.2, linewidth=0.4)
    if show_legend:
        ax.legend(fontsize=6, loc="upper left", framealpha=0.8)


# ---------------------------------------------------------------------------
# Panel B — DockQ speed
# ---------------------------------------------------------------------------

def draw_panel_b(ax: plt.Axes, rows: list[dict]) -> None:
    pairs = [r["pair"] for r in rows]
    pdb_cpp_times = [float(r["pdb_cpp_mean_s"]) for r in rows]
    dockq_times = [float(r["dockq_mean_s"]) for r in rows]
    speedups = [float(r["speedup_dockq_over_pdb_cpp"]) for r in rows]

    x = np.arange(len(pairs), dtype=float)
    width = 0.35

    ax.bar(
        x - width / 2,
        pdb_cpp_times,
        width=width,
        label="pdb_cpp",
        color=COLORS["pdb_cpp"],
        edgecolor="black",
        linewidth=0.3,
        alpha=0.9,
    )
    ax.bar(
        x + width / 2,
        dockq_times,
        width=width,
        label="DockQ",
        color=COLORS["DockQ"],
        edgecolor="black",
        linewidth=0.3,
        alpha=0.9,
    )

    # Annotate speedup above the DockQ bar
    for i, (dt, sp) in enumerate(zip(dockq_times, speedups)):
        ax.text(
            x[i] + width / 2,
            dt * 1.15,
            f"{sp:.1f}×",
            ha="center",
            va="bottom",
            fontsize=6,
            fontweight="bold",
        )

    ax.set_title("DockQ end-to-end", fontsize=8, fontweight="bold")
    ax.set_yscale("log")
    ax.set_ylabel("Execution time (s)", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [p.replace("_", "\n") for p in pairs],
        rotation=0,
        ha="center",
        fontsize=6,
    )
    ax.tick_params(axis="y", labelsize=6)
    ax.grid(True, axis="y", alpha=0.2, linewidth=0.4)
    ax.legend(fontsize=6, loc="upper left", framealpha=0.8)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _build_figure(common_data, dockq_rows, file_keys):
    """Create the multi-panel figure and return the Figure object."""
    n_common = len(file_keys)

    fig = plt.figure(figsize=(9.0, 5.5), constrained_layout=True)

    # Top row: one panel per protein size; bottom row: DockQ panel centered
    gs = fig.add_gridspec(2, n_common, height_ratios=[1, 0.85])
    top_axes = [fig.add_subplot(gs[0, i]) for i in range(n_common)]
    ax_dockq = fig.add_subplot(gs[1, n_common // 2])

    panel_letter = ord("A")
    for i, (ax, fk) in enumerate(zip(top_axes, file_keys)):
        label = SIZE_LABELS[i] if i < len(SIZE_LABELS) else f"File {i+1}"
        draw_common_panel(
            ax, common_data, fk,
            size_label=label,
            show_ylabel=(i == 0),
            show_legend=(i == 0),
        )
        ax.text(-0.08, 1.05, chr(panel_letter + i), transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="bottom")

    draw_panel_b(ax_dockq, dockq_rows)
    ax_dockq.text(-0.12, 1.05, chr(panel_letter + n_common),
                  transform=ax_dockq.transAxes,
                  fontsize=11, fontweight="bold", va="bottom")

    return fig


def main() -> None:
    args = parse_args()

    common_path = Path(args.common)
    dockq_path = Path(args.dockq)
    output_path = Path(args.output)

    if not common_path.exists():
        raise FileNotFoundError(f"Common CSV not found: {common_path}")
    if not dockq_path.exists():
        raise FileNotFoundError(f"DockQ CSV not found: {dockq_path}")

    common_data = load_common(common_path)
    dockq_rows = load_dockq(dockq_path)

    file_keys = _sorted_file_keys(common_data)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig = _build_figure(common_data, dockq_rows, file_keys)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved to: {output_path}")

    # Also save PNG copy for quick preview
    png_path = output_path.with_suffix(".png")
    fig2 = _build_figure(common_data, dockq_rows, file_keys)
    fig2.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig2)
    print(f"PNG preview saved to: {png_path}")


if __name__ == "__main__":
    main()
