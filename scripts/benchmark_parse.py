#!/usr/bin/env python3
"""Measure ordinary pdb_cpp coordinate-load throughput for input files."""

import argparse
import statistics
import time
from pathlib import Path

from pdb_cpp import Coor


def benchmark(path: Path, repeats: int) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    durations = []
    atom_count = 0
    for _ in range(repeats):
        start = time.perf_counter()
        coor = Coor(str(path))
        durations.append(time.perf_counter() - start)
        if coor.model_size() == 0:
            raise ValueError(f"No coordinate models loaded from {path}")
        atom_count = sum(model.size() for model in coor.get_all_Models())

    median_seconds = statistics.median(durations)
    atoms_per_second = atom_count / median_seconds if median_seconds else float("inf")
    print(
        f"{path}: {atom_count} atoms, median {median_seconds:.6f} s, "
        f"{atoms_per_second:.0f} atoms/s ({repeats} runs)"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path, help="PDB or mmCIF input files")
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be at least 1")
    for path in args.files:
        benchmark(path, args.repeats)


if __name__ == "__main__":
    main()