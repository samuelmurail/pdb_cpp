#!/usr/bin/env python3
# coding: utf-8

"""Benchmark coordinate read/write speed across structure libraries.

Targets:
- pdb_cpp
- pdb_numpy
- biopython (Bio.PDB)
- biotite

Usage examples
--------------
PYTHONPATH=src python benchmark/compare_io_speed.py
PYTHONPATH=src python benchmark/compare_io_speed.py --runs 10 --warmup 2
PYTHONPATH=src python benchmark/compare_io_speed.py --csv benchmark/io_speed.csv
"""

from __future__ import annotations

import argparse
import csv
import os
import time
from pathlib import Path
from statistics import mean, stdev
from typing import Any, Callable


Reader = Callable[[str], Any]
Writer = Callable[[Any, str], None]


DEFAULT_FILES = [
    "tests/input/1y0m.pdb",
    "tests/input/2rri.pdb",
    "tests/input/1y0m.cif",
    "tests/input/2rri.cif",
    "tests/input/9X0F.cif",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark read/write performance for coordinate files"
    )
    parser.add_argument("--runs", type=int, default=5, help="Measured runs per file/op")
    parser.add_argument("--warmup", type=int, default=1, help="Warmup runs per file/op")
    parser.add_argument(
        "--files",
        nargs="*",
        default=DEFAULT_FILES,
        help="Input coordinate files (PDB/mmCIF)",
    )
    parser.add_argument(
        "--csv",
        default="benchmark/io_speed_comparison.csv",
        help="Output CSV path",
    )
    return parser.parse_args()


def safe_import_backends() -> dict[str, dict[str, dict[str, Callable]]]:
    backends: dict[str, dict[str, dict[str, Callable]]] = {}

    try:
        from pdb_cpp import Coor

        def cpp_read(path: str):
            return Coor(path)

        def cpp_write(obj: Any, out_path: str):
            obj.write(out_path)

        backends["pdb_cpp"] = {
            "pdb": {"read": cpp_read, "write": cpp_write},
            "cif": {"read": cpp_read, "write": cpp_write},
        }
    except Exception:
        pass

    try:
        import pdb_numpy

        def numpy_read(path: str):
            return pdb_numpy.Coor(path)

        def numpy_write(obj: Any, out_path: str):
            obj.write(out_path, overwrite=True)

        backends["pdb_numpy"] = {
            "pdb": {"read": numpy_read, "write": numpy_write},
            "cif": {"read": numpy_read, "write": numpy_write},
        }
    except Exception:
        pass

    try:
        from Bio.PDB import MMCIFIO, MMCIFParser, PDBIO, PDBParser

        def biopython_read_pdb(path: str):
            parser = PDBParser(QUIET=True)
            return parser.get_structure("model", path)

        def biopython_write_pdb(obj: Any, out_path: str):
            io = PDBIO()
            io.set_structure(obj)
            io.save(out_path)

        def biopython_read_cif(path: str):
            parser = MMCIFParser(QUIET=True)
            return parser.get_structure("model", path)

        def biopython_write_cif(obj: Any, out_path: str):
            io = MMCIFIO()
            io.set_structure(obj)
            io.save(out_path)

        backends["biopython"] = {
            "pdb": {"read": biopython_read_pdb, "write": biopython_write_pdb},
            "cif": {"read": biopython_read_cif, "write": biopython_write_cif},
        }
    except Exception:
        pass

    try:
        import biotite.structure.io.pdb as pdb
        import biotite.structure.io.pdbx as pdbx

        def biotite_read_pdb(path: str):
            pdb_file = pdb.PDBFile.read(path)
            return pdb_file.get_structure(model=1)

        def biotite_write_pdb(obj: Any, out_path: str):
            pdb_file = pdb.PDBFile()
            pdb_file.set_structure(obj)
            pdb_file.write(out_path)

        def biotite_read_cif(path: str):
            cif_file = pdbx.CIFFile.read(path)
            return pdbx.get_structure(cif_file, model=1)

        def biotite_write_cif(obj: Any, out_path: str):
            cif_file = pdbx.CIFFile()
            pdbx.set_structure(cif_file, obj)
            cif_file.write(out_path)

        backends["biotite"] = {
            "pdb": {"read": biotite_read_pdb, "write": biotite_write_pdb},
            "cif": {"read": biotite_read_cif, "write": biotite_write_cif},
        }
    except Exception:
        pass

    return backends


def time_read(read_fn: Reader, path: str, warmup: int, runs: int) -> list[float]:
    for _ in range(warmup):
        read_fn(path)
    times = []
    for _ in range(runs):
        start = time.perf_counter()
        read_fn(path)
        times.append(time.perf_counter() - start)
    return times


def time_write(
    read_fn: Reader,
    write_fn: Writer,
    path: str,
    out_path: str,
    warmup: int,
    runs: int,
) -> list[float]:
    for _ in range(warmup):
        obj = read_fn(path)
        write_fn(obj, out_path)
    times = []
    for _ in range(runs):
        obj = read_fn(path)
        start = time.perf_counter()
        write_fn(obj, out_path)
        times.append(time.perf_counter() - start)
    return times


def summarize(times: list[float]) -> tuple[float, float]:
    if len(times) == 1:
        return times[0], 0.0
    return mean(times), stdev(times)


def main() -> None:
    args = parse_args()
    if args.runs < 1:
        raise ValueError("--runs must be >= 1")
    if args.warmup < 0:
        raise ValueError("--warmup must be >= 0")

    backends = safe_import_backends()
    if not backends:
        raise RuntimeError(
            "No backend import succeeded. Ensure environment has at least one of: "
            "pdb_cpp, pdb_numpy, biopython, biotite"
        )

    files = [Path(file_path) for file_path in args.files]
    for file_path in files:
        if not file_path.exists():
            raise FileNotFoundError(f"Input file not found: {file_path}")

    out_dir = Path("benchmark/tmp_io")
    out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []

    print(
        f"Benchmark settings: warmup={args.warmup}, runs={args.runs}, "
        f"files={len(files)}, backends={', '.join(sorted(backends))}"
    )
    print("library\tfile\toperation\tmean_s\tstd_s")

    for file_path in files:
        extension = file_path.suffix.lower().lstrip(".")
        if extension not in {"pdb", "cif"}:
            continue

        for backend_name, backend_impl in backends.items():
            if extension not in backend_impl:
                continue

            read_fn = backend_impl[extension]["read"]
            write_fn = backend_impl[extension]["write"]
            out_file = out_dir / f"{backend_name}_{file_path.stem}.{extension}"

            try:
                read_times = time_read(read_fn, str(file_path), args.warmup, args.runs)
                read_mean, read_std = summarize(read_times)
                print(
                    f"{backend_name}\t{file_path.name}\tread\t{read_mean:.6f}\t{read_std:.6f}"
                )
                rows.append(
                    {
                        "library": backend_name,
                        "file": str(file_path),
                        "format": extension,
                        "operation": "read",
                        "mean_s": read_mean,
                        "std_s": read_std,
                        "runs": args.runs,
                        "warmup": args.warmup,
                    }
                )

                write_times = time_write(
                    read_fn,
                    write_fn,
                    str(file_path),
                    str(out_file),
                    args.warmup,
                    args.runs,
                )
                write_mean, write_std = summarize(write_times)
                print(
                    f"{backend_name}\t{file_path.name}\twrite\t{write_mean:.6f}\t{write_std:.6f}"
                )
                rows.append(
                    {
                        "library": backend_name,
                        "file": str(file_path),
                        "format": extension,
                        "operation": "write",
                        "mean_s": write_mean,
                        "std_s": write_std,
                        "runs": args.runs,
                        "warmup": args.warmup,
                    }
                )
            except Exception as exc:
                print(f"{backend_name}\t{file_path.name}\tskipped\t{exc}")

            if out_file.exists():
                try:
                    os.remove(out_file)
                except OSError:
                    pass

    csv_path = Path(args.csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "library",
                "file",
                "format",
                "operation",
                "mean_s",
                "std_s",
                "runs",
                "warmup",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nCSV written to: {csv_path}")


if __name__ == "__main__":
    main()
