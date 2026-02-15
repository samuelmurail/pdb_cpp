#!/usr/bin/env python3
# coding: utf-8

"""Create grouped benchmark chart from CSV output.

Usage examples
--------------
python benchmark/plot_io_histogram.py
python benchmark/plot_io_histogram.py --input benchmark/io_speed_comparison.csv --output benchmark/io_speed_histogram.png
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot grouped operation chart with standard-error bars"
    )
    parser.add_argument(
        "--input",
        default="benchmark/io_speed_comparison.csv",
        help="Input CSV from compare_io_speed.py",
    )
    parser.add_argument(
        "--output",
        default="benchmark/io_speed_histogram.png",
        help="Output PNG path",
    )
    parser.add_argument(
        "--linear",
        action="store_true",
        help="Use linear y-axis instead of logarithmic",
    )
    return parser.parse_args()


def load_data(csv_path: Path) -> dict[str, dict[str, dict[str, list[float]]]]:
    values: dict[str, dict[str, dict[str, list[float]]]] = {}
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            library = row["library"]
            operation = row["operation"]
            file_path = row["file"]
            mean_s = float(row["mean_s"])
            values.setdefault(file_path, {}).setdefault(operation, {}).setdefault(
                library, []
            ).append(mean_s)
    return values


def mean_and_sem(values: list[float]) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    if len(values) == 1:
        return values[0], 0.0
    arr = np.asarray(values, dtype=float)
    mean_val = float(np.mean(arr))
    sem_val = float(np.std(arr, ddof=1) / math.sqrt(len(arr)))
    return mean_val, sem_val


def file_atom_count(file_path: Path) -> int:
    if not file_path.exists():
        return 0
    suffix = file_path.suffix.lower()
    count = 0
    with file_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if suffix == ".pdb":
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    count += 1
            elif suffix in {".cif", ".mmcif"}:
                stripped = line.lstrip()
                if stripped.startswith("ATOM ") or stripped.startswith("HETATM "):
                    count += 1
    return count


def main() -> None:
    args = parse_args()
    csv_path = Path(args.input)
    output_path = Path(args.output)

    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    data = load_data(csv_path)
    if not data:
        raise ValueError(f"No benchmark rows found in: {csv_path}")

    file_paths = [Path(file_str) for file_str in data.keys()]
    file_paths = [path for path in file_paths if path.exists()]
    if not file_paths:
        raise ValueError("No valid file paths found in CSV")

    file_paths.sort(key=lambda p: file_atom_count(p))

    operations = sorted(
        {
            operation
            for per_file in data.values()
            for operation in per_file.keys()
        }
    )
    libraries = sorted(
        {
            library
            for per_file in data.values()
            for per_op in per_file.values()
            for library in per_op.keys()
        }
    )

    preferred = ["pdb_cpp", "pdb_numpy", "biopython", "biotite"]
    libraries = [lib for lib in preferred if lib in libraries] + [
        lib for lib in libraries if lib not in preferred
    ]

    n_files = len(file_paths)
    fig, axes = plt.subplots(
        1, n_files, figsize=(5.2 * n_files, 5.8), constrained_layout=True
    )
    if n_files == 1:
        axes = [axes]

    size_labels = ["Small", "Medium", "Big"]

    for panel_index, (ax, file_path) in enumerate(zip(axes, file_paths)):
        file_key = str(file_path)
        x = np.arange(len(operations), dtype=float)
        width = 0.8 / max(1, len(libraries))

        for idx, library in enumerate(libraries):
            means = []
            sems = []
            for op in operations:
                op_values = data.get(file_key, {}).get(op, {}).get(library, [])
                m, se = mean_and_sem(op_values)
                means.append(m)
                sems.append(se)

            offsets = x - 0.4 + width / 2 + idx * width
            ax.bar(
                offsets,
                means,
                width=width,
                yerr=sems,
                capsize=3,
                label=library,
                alpha=0.9,
                edgecolor="black",
                linewidth=0.4,
            )

        pdb_id = file_path.stem.lower()
        atom_count = file_atom_count(file_path)
        if n_files == 3 and panel_index < 3:
            label = size_labels[panel_index]
            title = f"{label} – {pdb_id} ({atom_count} atoms)"
        else:
            title = f"{pdb_id} ({atom_count} atoms)"

        ax.set_xticks(x)
        ax.set_xticklabels(operations, rotation=25, ha="right")
        if not args.linear:
            ax.set_yscale("log")
        ax.set_title(title)
        ax.grid(True, axis="y", alpha=0.25)
        if panel_index == 0:
            ax.set_ylabel("Execution time (s)")
        ax.set_xlabel("Benchmark tests")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        borderaxespad=0.0,
        frameon=False,
    )
    fig.suptitle("Common benchmark by structure (mean ± standard error)")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Histogram written to: {output_path}")


if __name__ == "__main__":
    main()
