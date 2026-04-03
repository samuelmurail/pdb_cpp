#!/usr/bin/env python3
# coding: utf-8

"""Plot the results produced by compare_common_speed.py."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


COLORS = {
    "pdb_cpp": "#0072B2",
    "pdb_numpy": "#56B4E9",
    "biopython": "#E69F00",
    "biotite": "#009E73",
}
PREFERRED_LIB_ORDER = ["pdb_cpp", "pdb_numpy", "biopython", "biotite"]
OP_LABELS = {
    "read": "read",
    "write": "write",
    "select_within10_chainA": "select\n(within 10 A)",
    "get_aa_seq": "sequence",
    "rmsd_ca_shift": "RMSD",
    "dihedral_ca": "dihedral",
    "align_seq_chainA": "seq. align",
    "align_ca_self": "struct.\nalign",
    "hbond": "H-bonds",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot common benchmark results")
    parser.add_argument("--input", default="benchmark_common.csv", help="Input CSV file")
    parser.add_argument("--output", default="benchmark_common.png", help="Output plot path")
    return parser.parse_args()


def _atom_count(path: Path) -> int:
    if not path.exists():
        return 0
    suffix = path.suffix.lower()
    count = 0
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if suffix == ".pdb" and line.startswith(("ATOM", "HETATM")):
                count += 1
            elif suffix in {".cif", ".mmcif"} and line.lstrip().startswith(("ATOM ", "HETATM ")):
                count += 1
    return count


def _mean_sem(values: list[float]) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    if len(values) == 1:
        return values[0], 0.0
    arr = np.asarray(values, dtype=float)
    return float(np.mean(arr)), float(np.std(arr, ddof=1) / math.sqrt(len(arr)))


def load_common(csv_path: Path) -> dict[str, dict[str, dict[str, list[float]]]]:
    data: dict[str, dict[str, dict[str, list[float]]]] = {}
    with csv_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            file_key = row["file"]
            op = row["operation"]
            lib = row["library"]
            data.setdefault(file_key, {}).setdefault(op, {}).setdefault(lib, []).append(
                float(row["mean_s"])
            )
    return data


def draw_panel(ax: plt.Axes, data: dict, file_key: str, show_ylabel: bool) -> None:
    per_op = data[file_key]
    operations = list(per_op.keys())
    libraries = sorted(
        {lib for op_data in per_op.values() for lib in op_data},
        key=lambda lib: PREFERRED_LIB_ORDER.index(lib) if lib in PREFERRED_LIB_ORDER else 99,
    )

    x = np.arange(len(operations), dtype=float)
    width = 0.8 / max(1, len(libraries))
    for idx, lib in enumerate(libraries):
        means = []
        sems = []
        for op in operations:
            mean_value, sem_value = _mean_sem(per_op.get(op, {}).get(lib, []))
            means.append(mean_value)
            sems.append(sem_value)
        offsets = x - 0.4 + width / 2 + idx * width
        ax.bar(
            offsets,
            means,
            width=width,
            yerr=sems,
            capsize=2,
            label=lib,
            color=COLORS.get(lib, "#888888"),
            edgecolor="black",
            linewidth=0.3,
            alpha=0.9,
        )

    file_path = Path(file_key)
    atom_count = _atom_count(file_path)
    title = file_path.stem.upper()
    if atom_count:
        title += f" ({atom_count:,} atoms)"
    ax.set_title(title, fontsize=8, fontweight="bold")
    ax.set_yscale("log")
    if show_ylabel:
        ax.set_ylabel("Execution time (s)", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels([OP_LABELS.get(op, op) for op in operations], rotation=45, ha="right", fontsize=6)
    ax.tick_params(axis="y", labelsize=6)
    ax.grid(True, axis="y", alpha=0.2, linewidth=0.4)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    if not input_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    data = load_common(input_path)
    file_keys = sorted(data, key=lambda path: _atom_count(Path(path)))
    fig, axes = plt.subplots(1, len(file_keys), figsize=(4.2 * len(file_keys), 4.8), constrained_layout=True)
    if len(file_keys) == 1:
        axes = [axes]
    for idx, (ax, file_key) in enumerate(zip(axes, file_keys)):
        draw_panel(ax, data, file_key, show_ylabel=(idx == 0))
    axes[0].legend(fontsize=7, loc="upper left", framealpha=0.8)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot written to: {output_path}")


if __name__ == "__main__":
    main()