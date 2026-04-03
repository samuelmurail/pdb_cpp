#!/usr/bin/env python3
# coding: utf-8

"""Plot the results produced by compare_dockq_speed.py."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot DockQ benchmark results")
    parser.add_argument("--input", default="benchmark_dockq.csv", help="Input CSV file")
    parser.add_argument("--output", default="benchmark_dockq.png", help="Output plot path")
    return parser.parse_args()


def load_rows(csv_path: Path) -> list[dict[str, str]]:
    with csv_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    if not input_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    rows = load_rows(input_path)
    pairs = [row["pair"] for row in rows]
    dockq_times = [float(row["dockq_mean_s"]) for row in rows]
    pdb_cpp_times = [float(row["pdb_cpp_mean_s"]) for row in rows]
    speedups = [float(row["speedup_dockq_over_pdb_cpp"]) for row in rows]

    x = np.arange(len(pairs), dtype=float)
    width = 0.35
    fig, ax = plt.subplots(figsize=(max(6, 1.6 * len(pairs)), 4.8), constrained_layout=True)
    ax.bar(x - width / 2, pdb_cpp_times, width=width, label="pdb_cpp", color="#0072B2", edgecolor="black", linewidth=0.3)
    ax.bar(x + width / 2, dockq_times, width=width, label="DockQ", color="#D55E00", edgecolor="black", linewidth=0.3)

    for idx, (dockq_time, speedup) in enumerate(zip(dockq_times, speedups)):
        ax.text(x[idx] + width / 2, dockq_time * 1.12, f"{speedup:.1f}x", ha="center", va="bottom", fontsize=8, fontweight="bold")

    ax.set_title("DockQ benchmark", fontsize=10, fontweight="bold")
    ax.set_ylabel("Execution time (s)")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([pair.replace("_", "\n") for pair in pairs], fontsize=8)
    ax.grid(True, axis="y", alpha=0.2, linewidth=0.4)
    ax.legend(framealpha=0.8)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot written to: {output_path}")


if __name__ == "__main__":
    main()