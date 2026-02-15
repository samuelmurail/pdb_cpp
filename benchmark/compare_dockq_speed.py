#!/usr/bin/env python3
# coding: utf-8

"""Benchmark DockQ CLI vs pdb_cpp.analysis.dockQ.

Usage examples
--------------
python benchmark/compare_dockq_speed.py
python benchmark/compare_dockq_speed.py --runs 10 --warmup 2 --mode end-to-end
python benchmark/compare_dockq_speed.py --mode cached --csv benchmark/results_cached.csv
"""

from __future__ import annotations

import argparse
import csv
import statistics as stats
import subprocess
import time
from pathlib import Path
from typing import Any

from pdb_cpp import Coor, analysis


DEFAULT_PAIRS: list[dict[str, Any]] = [
    {
        "name": "1rxz_colabfold_vs_1rxz",
        "model": "tests/input/1rxz_colabfold_model_1.pdb",
        "native": "tests/input/1rxz.pdb",
        "dockq_args": [],
    },
    {
        "name": "model_vs_native",
        "model": "tests/input/model.pdb",
        "native": "tests/input/native.pdb",
        "dockq_args": ["--allowed_mismatches", "5"],
    },
    {
        "name": "1jd4_vs_5m6n",
        "model": "tests/input/1jd4.pdb",
        "native": "tests/input/5m6n.pdb",
        "dockq_args": ["--allowed_mismatches", "999", "--mapping", "AB:AB"],
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark DockQ CLI vs pdb_cpp DockQ")
    parser.add_argument("--runs", type=int, default=5, help="Measured repetitions per pair")
    parser.add_argument("--warmup", type=int, default=1, help="Warmup repetitions per pair")
    parser.add_argument(
        "--mode",
        choices=["end-to-end", "cached"],
        default="end-to-end",
        help=(
            "Benchmark mode: 'end-to-end' includes file parsing for pdb_cpp on every run; "
            "'cached' reuses preloaded Coor objects for pdb_cpp"
        ),
    )
    parser.add_argument(
        "--csv",
        default="benchmark/dockq_vs_pdb_cpp.csv",
        help="Output CSV path",
    )
    return parser.parse_args()


def run_dockq(model_path: str, native_path: str, dockq_args: list[str]) -> None:
    subprocess.run(
        ["DockQ", *dockq_args, model_path, native_path],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )


def run_pdb_cpp(
    model_path: str,
    native_path: str,
    mode: str,
    cached_model: Coor | None,
    cached_native: Coor | None,
) -> None:
    if mode == "cached":
        if cached_model is None or cached_native is None:
            raise ValueError("Cached mode requires preloaded coordinates")
        analysis.dockQ(cached_model, cached_native)
        return

    model_coor = Coor(model_path)
    native_coor = Coor(native_path)
    analysis.dockQ(model_coor, native_coor)


def benchmark_pair(pair: dict[str, Any], runs: int, warmup: int, mode: str) -> dict[str, Any]:
    model_path = pair["model"]
    native_path = pair["native"]
    dockq_args = pair["dockq_args"]

    cached_model = Coor(model_path) if mode == "cached" else None
    cached_native = Coor(native_path) if mode == "cached" else None

    for _ in range(warmup):
        run_dockq(model_path, native_path, dockq_args)
        run_pdb_cpp(model_path, native_path, mode, cached_model, cached_native)

    dockq_times = []
    cpp_times = []

    for _ in range(runs):
        t0 = time.perf_counter()
        run_dockq(model_path, native_path, dockq_args)
        dockq_times.append(time.perf_counter() - t0)

        t1 = time.perf_counter()
        run_pdb_cpp(model_path, native_path, mode, cached_model, cached_native)
        cpp_times.append(time.perf_counter() - t1)

    dockq_mean = stats.mean(dockq_times)
    dockq_median = stats.median(dockq_times)
    cpp_mean = stats.mean(cpp_times)
    cpp_median = stats.median(cpp_times)
    speedup = dockq_mean / cpp_mean if cpp_mean > 0 else float("inf")

    return {
        "pair": pair["name"],
        "dockq_mean_s": dockq_mean,
        "dockq_median_s": dockq_median,
        "pdb_cpp_mean_s": cpp_mean,
        "pdb_cpp_median_s": cpp_median,
        "speedup_dockq_over_pdb_cpp": speedup,
        "runs": runs,
        "warmup": warmup,
        "mode": mode,
    }


def write_csv(rows: list[dict[str, Any]], csv_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "pair",
        "dockq_mean_s",
        "dockq_median_s",
        "pdb_cpp_mean_s",
        "pdb_cpp_median_s",
        "speedup_dockq_over_pdb_cpp",
        "runs",
        "warmup",
        "mode",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def print_table(rows: list[dict[str, Any]]) -> None:
    print("pair\tDockQ_mean\tDockQ_median\tpdb_cpp_mean\tpdb_cpp_median\tspeedup")
    for row in rows:
        print(
            f"{row['pair']}\t"
            f"{row['dockq_mean_s']:.6f}\t"
            f"{row['dockq_median_s']:.6f}\t"
            f"{row['pdb_cpp_mean_s']:.6f}\t"
            f"{row['pdb_cpp_median_s']:.6f}\t"
            f"{row['speedup_dockq_over_pdb_cpp']:.2f}x"
        )


def main() -> None:
    args = parse_args()

    if args.runs < 1:
        raise ValueError("--runs must be >= 1")
    if args.warmup < 0:
        raise ValueError("--warmup must be >= 0")

    rows = [
        benchmark_pair(pair, runs=args.runs, warmup=args.warmup, mode=args.mode)
        for pair in DEFAULT_PAIRS
    ]

    print(
        f"Benchmark settings: warmup={args.warmup}, measured_runs={args.runs}, mode={args.mode}"
    )
    print_table(rows)

    csv_path = Path(args.csv)
    write_csv(rows, csv_path)
    print(f"\nCSV written to: {csv_path}")


if __name__ == "__main__":
    main()
