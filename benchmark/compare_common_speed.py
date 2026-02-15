#!/usr/bin/env python3
# coding: utf-8

"""Benchmark common structure operations across pdb_cpp/pdb_numpy/biopython/biotite.

Operations:
- read
- write
- select_within10_chainA
- get_aa_seq
- rmsd_ca_shift
- dihedral_ca
- align_seq_chainA
- align_ca_self

Usage
-----
PYTHONPATH=src python benchmark/compare_common_speed.py
PYTHONPATH=src python benchmark/compare_common_speed.py --runs 5 --warmup 1
PYTHONPATH=src python benchmark/compare_common_speed.py --files tests/input/1y0m.pdb tests/input/2rri.pdb
PYTHONPATH=src python benchmark/compare_common_speed.py --files tests/input/9X0F.cif
"""

from __future__ import annotations

import argparse
import csv
import os
import time
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, stdev
from typing import Any, Callable
import numpy as np


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
}


@dataclass
class Backend:
    name: str
    read: Callable[[str], Any]
    write: Callable[[Any, str], None]
    select_within10_chainA: Callable[[Any], int]
    get_aa_seq_len: Callable[[Any], int]
    get_ca_coords: Callable[[Any], np.ndarray]
    get_chain_seq_ca: Callable[[Any, str], tuple[str, np.ndarray]]
    align_ca_self: Callable[[str], float]


DEFAULT_FILES = [
    "tests/input/1y0m.pdb",
    "tests/input/1rxz.pdb",
    "tests/input/2mus.pdb",
    "tests/input/9X0F.cif",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark common structure operations for 4 Python packages"
    )
    parser.add_argument("--runs", type=int, default=5, help="Measured runs per operation")
    parser.add_argument("--warmup", type=int, default=1, help="Warmup runs per operation")
    parser.add_argument(
        "--files",
        nargs="*",
        default=DEFAULT_FILES,
        help="Input PDB files to benchmark",
    )
    parser.add_argument(
        "--csv",
        default="benchmark/common_speed_comparison.csv",
        help="Output CSV file",
    )
    return parser.parse_args()


def _safe_remove(path: str) -> None:
    try:
        if os.path.exists(path):
            os.remove(path)
    except OSError:
        pass


def _sum_seq_len_dict(seq_map: dict[str, str]) -> int:
    return sum(len(seq.replace("-", "")) for seq in seq_map.values())


def _sum_seq_len_list(seqs: list[str]) -> int:
    return sum(len(seq.replace("-", "")) for seq in seqs)


def _is_cif(path: str) -> bool:
    return Path(path).suffix.lower() in {".cif", ".mmcif"}


def _array_to_string(chars: Any) -> str:
    if isinstance(chars, str):
        return chars.strip()
    if isinstance(chars, (list, tuple)):
        return "".join(ch for ch in chars if ch not in ("\x00", " ")).strip()
    return str(chars).strip()


def _rmsd(coords_1: np.ndarray, coords_2: np.ndarray) -> float:
    n = min(coords_1.shape[0], coords_2.shape[0])
    if n == 0:
        return 0.0
    delta = coords_1[:n] - coords_2[:n]
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def _dihedral(points: np.ndarray) -> float:
    p0, p1, p2, p3 = points
    b0 = p1 - p0
    b1 = p2 - p1
    b2 = p3 - p2
    b1_norm = b1 / (np.linalg.norm(b1) + 1e-12)
    v = b0 - np.dot(b0, b1_norm) * b1_norm
    w = b2 - np.dot(b2, b1_norm) * b1_norm
    x = np.dot(v, w)
    y = np.dot(np.cross(b1_norm, v), w)
    return float(np.arctan2(y, x))


def _kabsch_rmsd(coords_1: np.ndarray, coords_2: np.ndarray) -> float:
    n = min(coords_1.shape[0], coords_2.shape[0])
    if n == 0:
        return 0.0
    p = coords_1[:n].astype(float)
    q = coords_2[:n].astype(float)
    p_cent = p.mean(axis=0)
    q_cent = q.mean(axis=0)
    p0 = p - p_cent
    q0 = q - q_cent
    h = p0.T @ q0
    u, _, vt = np.linalg.svd(h)
    if np.linalg.det(vt.T @ u.T) < 0:
        vt[-1, :] *= -1
    r = vt.T @ u.T
    p_aligned = p0 @ r
    return _rmsd(p_aligned, q0)


def _sequence_alignment_pairs(seq_1: str, seq_2: str) -> list[tuple[int, int]]:
    try:
        from Bio import pairwise2

        alignment = pairwise2.align.globalms(
            seq_1, seq_2, 2, -1, -1, -0.5, one_alignment_only=True
        )[0]
        a1, a2 = alignment.seqA, alignment.seqB
        pairs: list[tuple[int, int]] = []
        i1 = 0
        i2 = 0
        for c1, c2 in zip(a1, a2):
            if c1 != "-" and c2 != "-":
                pairs.append((i1, i2))
            if c1 != "-":
                i1 += 1
            if c2 != "-":
                i2 += 1
        return pairs
    except Exception:
        n = min(len(seq_1), len(seq_2))
        return [(i, i) for i in range(n)]


def _align_seq_rmsd(seq_1: str, coords_1: np.ndarray, seq_2: str, coords_2: np.ndarray) -> float:
    pairs = _sequence_alignment_pairs(seq_1, seq_2)
    if not pairs:
        return 0.0
    p = np.asarray([coords_1[i] for i, _ in pairs], dtype=float)
    q = np.asarray([coords_2[j] for _, j in pairs], dtype=float)
    return _kabsch_rmsd(p, q)


def build_backends() -> list[Backend]:
    backends: list[Backend] = []

    try:
        from pdb_cpp import Coor, core

        def cpp_read(path: str):
            return Coor(path)

        def cpp_write(obj: Any, out_path: str):
            obj.write(out_path)

        def cpp_select_within10_chainA(obj: Any) -> int:
            return obj.select_atoms("within 10.0 of chain A").size()

        def cpp_get_aa_seq_len(obj: Any) -> int:
            if hasattr(obj, "get_aa_seq"):
                return _sum_seq_len_dict(obj.get_aa_seq())
            if hasattr(obj, "get_aa_sequences"):
                return _sum_seq_len_list(obj.get_aa_sequences())
            raise AttributeError("No sequence API found on pdb_cpp Coor")

        def cpp_get_ca_coords(obj: Any) -> np.ndarray:
            model = obj.get_Models(0)
            idx = obj.get_index_select("name CA")
            if len(idx) == 0:
                return np.zeros((0, 3), dtype=float)
            sel = np.asarray(idx, dtype=int)
            xs = np.asarray(model.get_x(), dtype=float)
            ys = np.asarray(model.get_y(), dtype=float)
            zs = np.asarray(model.get_z(), dtype=float)
            return np.stack([xs[sel], ys[sel], zs[sel]], axis=1)

        def cpp_get_chain_seq_ca(obj: Any, chain_id: str) -> tuple[str, np.ndarray]:
            model = obj.get_Models(0)
            names = model.get_name()
            chains = model.get_chain()
            resnames = model.get_resname()
            xs = model.get_x()
            ys = model.get_y()
            zs = model.get_z()
            seq = []
            coords = []
            for i in range(len(names)):
                if _array_to_string(names[i]) != "CA":
                    continue
                if _array_to_string(chains[i]) != chain_id:
                    continue
                seq.append(AA3_TO_1.get(_array_to_string(resnames[i]), "X"))
                coords.append([float(xs[i]), float(ys[i]), float(zs[i])])
            if not coords:
                return "", np.zeros((0, 3), dtype=float)
            return "".join(seq), np.asarray(coords, dtype=float)

        def cpp_align_ca_self(path: str) -> float:
            model = Coor(path)
            native = Coor(path)
            idx_model = model.get_index_select("name CA")
            idx_native = native.get_index_select("name CA")
            if len(idx_model) == 0 or len(idx_native) == 0:
                return 0.0
            n = min(len(idx_model), len(idx_native))
            idx_model = idx_model[:n]
            idx_native = idx_native[:n]
            return core.coor_align(model, native, idx_model, idx_native)

        backends.append(
            Backend(
                name="pdb_cpp",
                read=cpp_read,
                write=cpp_write,
                select_within10_chainA=cpp_select_within10_chainA,
                get_aa_seq_len=cpp_get_aa_seq_len,
                get_ca_coords=cpp_get_ca_coords,
                get_chain_seq_ca=cpp_get_chain_seq_ca,
                align_ca_self=cpp_align_ca_self,
            )
        )
    except Exception:
        pass

    try:
        import pdb_numpy
        from pdb_numpy import alignement

        def numpy_read(path: str):
            return pdb_numpy.Coor(path)

        def numpy_write(obj: Any, out_path: str):
            obj.write(out_path, overwrite=True)

        def numpy_select_within10_chainA(obj: Any) -> int:
            return obj.select_atoms("within 10.0 of chain A").len

        def numpy_get_aa_seq_len(obj: Any) -> int:
            if hasattr(obj, "get_aa_seq"):
                return _sum_seq_len_dict(obj.get_aa_seq())
            if hasattr(obj, "get_aa_sequences"):
                return _sum_seq_len_list(obj.get_aa_sequences())
            raise AttributeError("No sequence API found on pdb_numpy Coor")

        def numpy_get_ca_coords(obj: Any) -> np.ndarray:
            idx = obj.get_index_select("name CA")
            if len(idx) == 0:
                return np.zeros((0, 3), dtype=float)
            model = obj.models[obj.active_model]
            sel = np.asarray(idx, dtype=int)
            return np.stack(
                [
                    np.asarray(model.x, dtype=float)[sel],
                    np.asarray(model.y, dtype=float)[sel],
                    np.asarray(model.z, dtype=float)[sel],
                ],
                axis=1,
            )

        def numpy_get_chain_seq_ca(obj: Any, chain_id: str) -> tuple[str, np.ndarray]:
            model = obj.models[obj.active_model]
            names = np.asarray(model.name)
            chains = np.asarray(model.chain)
            resnames = np.asarray(model.resname)
            mask = (names == "CA") & (chains == chain_id)
            if not np.any(mask):
                return "", np.zeros((0, 3), dtype=float)
            seq = "".join(AA3_TO_1.get(str(r), "X") for r in resnames[mask])
            coords = np.stack(
                [
                    np.asarray(model.x, dtype=float)[mask],
                    np.asarray(model.y, dtype=float)[mask],
                    np.asarray(model.z, dtype=float)[mask],
                ],
                axis=1,
            )
            return seq, coords

        def numpy_align_ca_self(path: str) -> float:
            model = pdb_numpy.Coor(path)
            native = pdb_numpy.Coor(path)
            idx_model = model.get_index_select("name CA")
            idx_native = native.get_index_select("name CA")
            if len(idx_model) == 0 or len(idx_native) == 0:
                return 0.0
            n = min(len(idx_model), len(idx_native))
            idx_model = idx_model[:n]
            idx_native = idx_native[:n]
            return alignement.coor_align(model, native, idx_model, idx_native)

        backends.append(
            Backend(
                name="pdb_numpy",
                read=numpy_read,
                write=numpy_write,
                select_within10_chainA=numpy_select_within10_chainA,
                get_aa_seq_len=numpy_get_aa_seq_len,
                get_ca_coords=numpy_get_ca_coords,
                get_chain_seq_ca=numpy_get_chain_seq_ca,
                align_ca_self=numpy_align_ca_self,
            )
        )
    except Exception:
        pass

    try:
        import numpy as np
        from Bio.PDB import MMCIFIO, MMCIFParser, PDBIO, PDBParser, Superimposer
        from Bio.PDB.Polypeptide import is_aa

        def bio_read(path: str):
            parser = MMCIFParser(QUIET=True) if _is_cif(path) else PDBParser(QUIET=True)
            return parser.get_structure("model", path)

        def bio_write(obj: Any, out_path: str):
            io = MMCIFIO() if _is_cif(out_path) else PDBIO()
            io.set_structure(obj)
            io.save(out_path)

        def bio_select_within10_chainA(obj: Any) -> int:
            chain_a_coords = [
                atom.coord
                for atom in obj.get_atoms()
                if atom.get_parent().get_parent().id == "A"
            ]
            if not chain_a_coords:
                return 0
            ref = np.asarray(chain_a_coords, dtype=float)
            all_coords = [atom.coord for atom in obj.get_atoms()]
            all_arr = np.asarray(all_coords, dtype=float)

            threshold2 = 100.0
            selected = 0
            chunk = 256
            for start in range(0, all_arr.shape[0], chunk):
                block = all_arr[start : start + chunk]
                d2 = np.sum((block[:, None, :] - ref[None, :, :]) ** 2, axis=2)
                selected += int(np.count_nonzero(np.any(d2 <= threshold2, axis=1)))
            return selected

        def bio_get_aa_seq_len(obj: Any) -> int:
            total = 0
            for model in obj:
                for chain in model:
                    for residue in chain:
                        if is_aa(residue, standard=False):
                            total += 1
            return total

        def bio_get_ca_coords(obj: Any) -> np.ndarray:
            coords = [atom.coord for atom in obj.get_atoms() if atom.get_name() == "CA"]
            if not coords:
                return np.zeros((0, 3), dtype=float)
            return np.asarray(coords, dtype=float)

        def bio_get_chain_seq_ca(obj: Any, chain_id: str) -> tuple[str, np.ndarray]:
            seq = []
            coords = []
            for atom in obj.get_atoms():
                if atom.get_name() != "CA":
                    continue
                residue = atom.get_parent()
                chain = residue.get_parent()
                if chain.id != chain_id:
                    continue
                seq.append(AA3_TO_1.get(residue.get_resname(), "X"))
                coords.append(atom.coord)
            if not coords:
                return "", np.zeros((0, 3), dtype=float)
            return "".join(seq), np.asarray(coords, dtype=float)

        def bio_align_ca_self(path: str) -> float:
            fixed = bio_read(path)
            mobile = bio_read(path)
            fixed_ca = [atom for atom in fixed.get_atoms() if atom.get_name() == "CA"]
            mobile_ca = [atom for atom in mobile.get_atoms() if atom.get_name() == "CA"]
            n = min(len(fixed_ca), len(mobile_ca))
            if n == 0:
                return 0.0
            sup = Superimposer()
            sup.set_atoms(fixed_ca[:n], mobile_ca[:n])
            return float(sup.rms)

        backends.append(
            Backend(
                name="biopython",
                read=bio_read,
                write=bio_write,
                select_within10_chainA=bio_select_within10_chainA,
                get_aa_seq_len=bio_get_aa_seq_len,
                get_ca_coords=bio_get_ca_coords,
                get_chain_seq_ca=bio_get_chain_seq_ca,
                align_ca_self=bio_align_ca_self,
            )
        )
    except Exception:
        pass

    try:
        import numpy as np
        import biotite.structure as struc
        import biotite.structure.io.pdb as pdb
        import biotite.structure.io.pdbx as pdbx

        def bt_read(path: str):
            if _is_cif(path):
                cif_file = pdbx.CIFFile.read(path)
                return pdbx.get_structure(cif_file, model=1)
            pdb_file = pdb.PDBFile.read(path)
            return pdb_file.get_structure(model=1)

        def bt_write(obj: Any, out_path: str):
            if _is_cif(out_path):
                cif_file = pdbx.CIFFile()
                pdbx.set_structure(cif_file, obj)
                cif_file.write(out_path)
                return
            pdb_file = pdb.PDBFile()
            pdb_file.set_structure(obj)
            pdb_file.write(out_path)

        def bt_select_within10_chainA(obj: Any) -> int:
            ref = obj[obj.chain_id == "A"]
            if ref.array_length() == 0:
                return 0
            all_coords = obj.coord
            ref_coords = ref.coord
            threshold2 = 100.0
            selected = 0
            chunk = 256
            for start in range(0, all_coords.shape[0], chunk):
                block = all_coords[start : start + chunk]
                d2 = np.sum((block[:, None, :] - ref_coords[None, :, :]) ** 2, axis=2)
                selected += int(np.count_nonzero(np.any(d2 <= threshold2, axis=1)))
            return selected

        def bt_get_aa_seq_len(obj: Any) -> int:
            _, res_names = struc.get_residues(obj)
            return sum(1 for name in res_names if name in AA3_TO_1)

        def bt_get_ca_coords(obj: Any) -> np.ndarray:
            mask = obj.atom_name == "CA"
            if not np.any(mask):
                return np.zeros((0, 3), dtype=float)
            return np.asarray(obj.coord[mask], dtype=float)

        def bt_get_chain_seq_ca(obj: Any, chain_id: str) -> tuple[str, np.ndarray]:
            mask = (obj.atom_name == "CA") & (obj.chain_id == chain_id)
            if not np.any(mask):
                return "", np.zeros((0, 3), dtype=float)
            seq = "".join(AA3_TO_1.get(str(r), "X") for r in obj.res_name[mask])
            return seq, np.asarray(obj.coord[mask], dtype=float)

        def bt_align_ca_self(path: str) -> float:
            fixed = bt_read(path)
            mobile = bt_read(path)
            fixed_ca = fixed[fixed.atom_name == "CA"]
            mobile_ca = mobile[mobile.atom_name == "CA"]
            n = min(fixed_ca.array_length(), mobile_ca.array_length())
            if n == 0:
                return 0.0
            fitted, _ = struc.superimpose(fixed_ca[:n], mobile_ca[:n])
            return float(struc.rmsd(fixed_ca[:n], fitted))

        backends.append(
            Backend(
                name="biotite",
                read=bt_read,
                write=bt_write,
                select_within10_chainA=bt_select_within10_chainA,
                get_aa_seq_len=bt_get_aa_seq_len,
                get_ca_coords=bt_get_ca_coords,
                get_chain_seq_ca=bt_get_chain_seq_ca,
                align_ca_self=bt_align_ca_self,
            )
        )
    except Exception:
        pass

    return backends


def benchmark_call(fn: Callable[[], Any], warmup: int, runs: int) -> tuple[float, float]:
    for _ in range(warmup):
        fn()

    times = []
    for _ in range(runs):
        start = time.perf_counter()
        fn()
        times.append(time.perf_counter() - start)

    if len(times) == 1:
        return times[0], 0.0
    return mean(times), stdev(times)


def main() -> None:
    args = parse_args()
    if args.runs < 1:
        raise ValueError("--runs must be >= 1")
    if args.warmup < 0:
        raise ValueError("--warmup must be >= 0")

    files = [Path(file_path) for file_path in args.files]
    for file_path in files:
        if not file_path.exists():
            raise FileNotFoundError(f"Missing file: {file_path}")
        if file_path.suffix.lower() not in {".pdb", ".cif", ".mmcif"}:
            raise ValueError(f"Only PDB/mmCIF files are supported for common-op benchmark: {file_path}")

    backends = build_backends()
    if not backends:
        raise RuntimeError("No backend could be imported")

    backend_names = [backend.name for backend in backends]
    print(
        f"Benchmark settings: warmup={args.warmup}, runs={args.runs}, files={len(files)}, "
        f"backends={', '.join(backend_names)}"
    )

    operations = [
        "read",
        "write",
        "select_within10_chainA",
        "get_aa_seq",
        "rmsd_ca_shift",
        "dihedral_ca",
        "align_seq_chainA",
        "align_ca_self",
    ]
    rows: list[dict[str, Any]] = []

    print("library\tfile\toperation\tmean_s\tstd_s\tspeedup_vs_pdb_cpp\tpdb_cpp_superior")

    for file_path in files:
        per_file_op_times: dict[tuple[str, str], float] = {}

        for backend in backends:
            out_path = f"benchmark/tmp_io/{backend.name}_{file_path.name}"
            backend_key_list: list[tuple[str, str]] = []
            backend_row_indices: list[int] = []
            try:
                read_mean, read_std = benchmark_call(
                    lambda b=backend, p=str(file_path): b.read(p), args.warmup, args.runs
                )
                per_file_op_times[(backend.name, "read")] = read_mean
                backend_key_list.append((backend.name, "read"))
                rows.append(
                    {
                        "library": backend.name,
                        "file": str(file_path),
                        "operation": "read",
                        "mean_s": read_mean,
                        "std_s": read_std,
                        "runs": args.runs,
                        "warmup": args.warmup,
                    }
                )
                backend_row_indices.append(len(rows) - 1)

                write_mean, write_std = benchmark_call(
                    lambda b=backend, p=str(file_path), o=out_path: b.write(b.read(p), o),
                    args.warmup,
                    args.runs,
                )
                per_file_op_times[(backend.name, "write")] = write_mean
                backend_key_list.append((backend.name, "write"))
                rows.append(
                    {
                        "library": backend.name,
                        "file": str(file_path),
                        "operation": "write",
                        "mean_s": write_mean,
                        "std_s": write_std,
                        "runs": args.runs,
                        "warmup": args.warmup,
                    }
                )
                backend_row_indices.append(len(rows) - 1)

                obj = backend.read(str(file_path))
                select_mean, select_std = benchmark_call(
                    lambda b=backend, current=obj: b.select_within10_chainA(current),
                    args.warmup,
                    args.runs,
                )
                per_file_op_times[(backend.name, "select_within10_chainA")] = select_mean
                backend_key_list.append((backend.name, "select_within10_chainA"))
                rows.append(
                    {
                        "library": backend.name,
                        "file": str(file_path),
                        "operation": "select_within10_chainA",
                        "mean_s": select_mean,
                        "std_s": select_std,
                        "runs": args.runs,
                        "warmup": args.warmup,
                    }
                )
                backend_row_indices.append(len(rows) - 1)

                seq_mean, seq_std = benchmark_call(
                    lambda b=backend, current=obj: b.get_aa_seq_len(current),
                    args.warmup,
                    args.runs,
                )
                per_file_op_times[(backend.name, "get_aa_seq")] = seq_mean
                backend_key_list.append((backend.name, "get_aa_seq"))
                rows.append(
                    {
                        "library": backend.name,
                        "file": str(file_path),
                        "operation": "get_aa_seq",
                        "mean_s": seq_mean,
                        "std_s": seq_std,
                        "runs": args.runs,
                        "warmup": args.warmup,
                    }
                )
                backend_row_indices.append(len(rows) - 1)

                rmsd_mean, rmsd_std = benchmark_call(
                    lambda b=backend, current=obj: _rmsd(
                        b.get_ca_coords(current),
                        b.get_ca_coords(current) + np.array([1.0, 0.0, 0.0]),
                    ),
                    args.warmup,
                    args.runs,
                )
                per_file_op_times[(backend.name, "rmsd_ca_shift")] = rmsd_mean
                backend_key_list.append((backend.name, "rmsd_ca_shift"))
                rows.append(
                    {
                        "library": backend.name,
                        "file": str(file_path),
                        "operation": "rmsd_ca_shift",
                        "mean_s": rmsd_mean,
                        "std_s": rmsd_std,
                        "runs": args.runs,
                        "warmup": args.warmup,
                    }
                )
                backend_row_indices.append(len(rows) - 1)

                dihedral_mean, dihedral_std = benchmark_call(
                    lambda b=backend, current=obj: (
                        _dihedral(b.get_ca_coords(current)[:4])
                        if b.get_ca_coords(current).shape[0] >= 4
                        else 0.0
                    ),
                    args.warmup,
                    args.runs,
                )
                per_file_op_times[(backend.name, "dihedral_ca")] = dihedral_mean
                backend_key_list.append((backend.name, "dihedral_ca"))
                rows.append(
                    {
                        "library": backend.name,
                        "file": str(file_path),
                        "operation": "dihedral_ca",
                        "mean_s": dihedral_mean,
                        "std_s": dihedral_std,
                        "runs": args.runs,
                        "warmup": args.warmup,
                    }
                )
                backend_row_indices.append(len(rows) - 1)

                align_seq_mean, align_seq_std = benchmark_call(
                    lambda b=backend, current=obj: (
                        lambda seq_1, coords_1, seq_2, coords_2: _align_seq_rmsd(
                            seq_1, coords_1, seq_2, coords_2
                        )
                    )(*b.get_chain_seq_ca(current, "A"), *b.get_chain_seq_ca(current, "A")),
                    args.warmup,
                    args.runs,
                )
                per_file_op_times[(backend.name, "align_seq_chainA")] = align_seq_mean
                backend_key_list.append((backend.name, "align_seq_chainA"))
                rows.append(
                    {
                        "library": backend.name,
                        "file": str(file_path),
                        "operation": "align_seq_chainA",
                        "mean_s": align_seq_mean,
                        "std_s": align_seq_std,
                        "runs": args.runs,
                        "warmup": args.warmup,
                    }
                )
                backend_row_indices.append(len(rows) - 1)

                align_mean, align_std = benchmark_call(
                    lambda b=backend, p=str(file_path): b.align_ca_self(p),
                    args.warmup,
                    args.runs,
                )
                per_file_op_times[(backend.name, "align_ca_self")] = align_mean
                backend_key_list.append((backend.name, "align_ca_self"))
                rows.append(
                    {
                        "library": backend.name,
                        "file": str(file_path),
                        "operation": "align_ca_self",
                        "mean_s": align_mean,
                        "std_s": align_std,
                        "runs": args.runs,
                        "warmup": args.warmup,
                    }
                )
                backend_row_indices.append(len(rows) - 1)
            except Exception as exc:
                for key in backend_key_list:
                    per_file_op_times.pop(key, None)
                for index in sorted(backend_row_indices, reverse=True):
                    if 0 <= index < len(rows):
                        rows.pop(index)
                print(f"{backend.name}\t{file_path.name}\tskipped\tNA\tNA\tNA\tNA ({exc})")
            finally:
                _safe_remove(out_path)

        for backend in backends:
            for op in operations:
                key = (backend.name, op)
                if key not in per_file_op_times:
                    continue
                current = per_file_op_times[key]
                base = per_file_op_times.get(("pdb_cpp", op))
                speedup = base / current if base and current > 0 else None
                superior = "NA"
                if base is not None:
                    superior = "YES" if base <= current else "NO"
                speedup_txt = "NA" if speedup is None else f"{speedup:.2f}x"
                matching = [
                    item["std_s"]
                    for item in rows
                    if item["library"] == backend.name
                    and item["file"] == str(file_path)
                    and item["operation"] == op
                ]
                if not matching:
                    continue
                std = matching[0]
                print(
                    f"{backend.name}\t{file_path.name}\t{op}\t{current:.6f}\t{std:.6f}\t{speedup_txt}\t{superior}"
                )

    csv_path = Path(args.csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "library",
                "file",
                "operation",
                "mean_s",
                "std_s",
                "runs",
                "warmup",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print("\nSummary against pdb_cpp (aggregated over files):")
    for op in operations:
        op_rows = [row for row in rows if row["operation"] == op]
        if not op_rows:
            print(f"- {op}: no data")
            continue
        cpp_means = [row["mean_s"] for row in op_rows if row["library"] == "pdb_cpp"]
        if cpp_means:
            cpp_mean = mean(cpp_means)
            op_means_by_lib = {}
            for row in op_rows:
                op_means_by_lib.setdefault(row["library"], []).append(row["mean_s"])
            fastest = min((mean(vals), lib) for lib, vals in op_means_by_lib.items())[1]
            verdict = "pdb_cpp fastest" if fastest == "pdb_cpp" else f"{fastest} fastest"
            print(f"- {op}: {verdict} (pdb_cpp avg={cpp_mean:.6f}s)")
        else:
            print(f"- {op}: pdb_cpp not available")

    print(f"\nCSV written to: {csv_path}")


if __name__ == "__main__":
    main()
