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

The script benchmarks a fixed set of operations on a user-provided list of
structures and writes one CSV row per backend, file, and operation.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
import time
import zlib
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
    rmsd_ca_shift: Callable[[Any, Any], float]
    dihedral_ca: Callable[[Any], float]
    align_seq_chainA: Callable[[str, str], float]
    align_ca_pair: Callable[[str, str], float]
    hbond_count: Callable[[Any], int]


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
        default="benchmark_common.csv",
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


def _report_backend_skip(name: str, exc: Exception) -> None:
    print(f"Skipping backend '{name}': {exc}", file=sys.stderr)


def _random_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    quat = rng.normal(size=4)
    quat /= np.linalg.norm(quat)
    w, x, y, z = quat
    return np.asarray(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ],
        dtype=float,
    )


def _create_transformed_structure(input_path: Path, output_path: Path) -> None:
    from Bio.PDB import MMCIFIO, MMCIFParser, PDBIO, PDBParser

    rng = np.random.default_rng(zlib.crc32(str(input_path).encode("utf-8")))
    rotation = _random_rotation_matrix(rng)
    translation = rng.uniform(-25.0, 25.0, size=3)

    parser = MMCIFParser(QUIET=True) if _is_cif(str(input_path)) else PDBParser(QUIET=True)
    structure = parser.get_structure("transformed", str(input_path))
    for atom in structure.get_atoms():
        atom.transform(rotation, translation)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    io = MMCIFIO() if _is_cif(str(output_path)) else PDBIO()
    io.set_structure(structure)
    io.save(str(output_path))


def _dihedral_numpy(pts: np.ndarray) -> float:
    """Pure-NumPy fallback: one dihedral from the first 4 rows of pts."""
    p0, p1, p2, p3 = pts[0], pts[1], pts[2], pts[3]
    b0 = p1 - p0; b1 = p2 - p1; b2 = p3 - p2
    b1n = b1 / (np.linalg.norm(b1) + 1e-12)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def _dihedral(pts: np.ndarray) -> float:
    """Dihedral angle (degrees) from 4 points; delegates to C++ when available."""
    try:
        from pdb_cpp.core import compute_dihedrals as _cpp_dihed
        arr = np.asarray(pts[:4], dtype=np.float32)
        result = _cpp_dihed(arr)
        return float(result[0]) if result.shape[0] > 0 else 0.0
    except Exception:
        return _dihedral_numpy(pts)


def build_backends() -> list[Backend]:
    backends: list[Backend] = []

    try:
        from pdb_cpp import Coor, analysis as pdb_analysis, core, geom as pdb_geom
        from pdb_cpp import hbond as pdb_hbond

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

        def cpp_rmsd_ca_shift(reference_obj: Any, mobile_obj: Any) -> float:
            return float(pdb_analysis.rmsd(mobile_obj, reference_obj, selection="name CA")[0])

        def cpp_dihedral_ca(obj: Any) -> float:
            ca = obj.select_atoms("name CA")
            if ca.len < 4:
                return 0.0
            return float(pdb_geom.compute_dihedrals(np.asarray(ca.xyz[:4], dtype=np.float32))[0])

        def cpp_align_seq_chainA(reference_path: str, mobile_path: str) -> float:
            native = Coor(reference_path)
            model = Coor(mobile_path)
            rmsds, _, _ = core.align_seq_based(
                model,
                native,
                chain_1=["A"],
                chain_2=["A"],
                back_names=["CA"],
            )
            return float(rmsds[0]) if rmsds else 0.0

        def cpp_align_ca_pair(reference_path: str, mobile_path: str) -> float:
            native = Coor(reference_path)
            model = Coor(mobile_path)
            idx_model = model.get_index_select("name CA")
            idx_native = native.get_index_select("name CA")
            if len(idx_model) == 0 or len(idx_native) == 0:
                return 0.0
            n = min(len(idx_model), len(idx_native))
            idx_model = idx_model[:n]
            idx_native = idx_native[:n]
            return core.coor_align(model, native, idx_model, idx_native)

        def cpp_hbond_count(obj: Any) -> int:
            return len(pdb_hbond.hbonds(obj, angle_cutoff=120.0)[0])

        backends.append(
            Backend(
                name="pdb_cpp",
                read=cpp_read,
                write=cpp_write,
                select_within10_chainA=cpp_select_within10_chainA,
                get_aa_seq_len=cpp_get_aa_seq_len,
                rmsd_ca_shift=cpp_rmsd_ca_shift,
                dihedral_ca=cpp_dihedral_ca,
                align_seq_chainA=cpp_align_seq_chainA,
                align_ca_pair=cpp_align_ca_pair,
                hbond_count=cpp_hbond_count,
            )
        )
    except Exception as exc:
        _report_backend_skip("pdb_cpp", exc)
        pass

    try:
        import pdb_numpy
        import pdb_numpy.analysis as numpy_analysis
        from pdb_numpy import alignement, DSSP as pdb_numpy_DSSP

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

        def numpy_rmsd_ca_shift(reference_obj: Any, mobile_obj: Any) -> float:
            rmsds = numpy_analysis.rmsd(mobile_obj, reference_obj, selection="name CA")
            return float(rmsds[0]) if rmsds else 0.0

        def numpy_dihedral_ca(obj: Any) -> float:
            ca = obj.select_atoms("name CA")
            if ca.len < 4:
                return 0.0
            model = ca.models[ca.active_model]
            pts = np.asarray(
                [
                    [float(model.x[i]), float(model.y[i]), float(model.z[i])]
                    for i in range(min(4, model.len))
                ],
                dtype=float,
            )
            return _dihedral_numpy(pts) if pts.shape[0] >= 4 else 0.0

        def numpy_align_seq_chainA(reference_path: str, mobile_path: str) -> float:
            native = pdb_numpy.Coor(reference_path)
            model = pdb_numpy.Coor(mobile_path)
            rmsds, _ = alignement.rmsd_seq_based(
                model,
                native,
                chain_1=["A"],
                chain_2=["A"],
                back_names=["CA"],
            )
            return float(rmsds[0]) if rmsds else 0.0

        def numpy_align_ca_pair(reference_path: str, mobile_path: str) -> float:
            native = pdb_numpy.Coor(reference_path)
            model = pdb_numpy.Coor(mobile_path)
            idx_model = model.get_index_select("name CA")
            idx_native = native.get_index_select("name CA")
            if len(idx_model) == 0 or len(idx_native) == 0:
                return 0.0
            n = min(len(idx_model), len(idx_native))
            idx_model = idx_model[:n]
            idx_native = idx_native[:n]
            return alignement.coor_align(model, native, idx_model, idx_native)

        def numpy_hbond_count(obj: Any) -> int:
            mat = pdb_numpy_DSSP.compute_Hbond_matrix(obj.models[obj.active_model])
            return int(np.sum(mat))

        backends.append(
            Backend(
                name="pdb_numpy",
                read=numpy_read,
                write=numpy_write,
                select_within10_chainA=numpy_select_within10_chainA,
                get_aa_seq_len=numpy_get_aa_seq_len,
                rmsd_ca_shift=numpy_rmsd_ca_shift,
                dihedral_ca=numpy_dihedral_ca,
                align_seq_chainA=numpy_align_seq_chainA,
                align_ca_pair=numpy_align_ca_pair,
                hbond_count=numpy_hbond_count,
            )
        )
    except Exception as exc:
        _report_backend_skip("pdb_numpy", exc)
        pass

    try:
        from Bio import Align
        from Bio.PDB import MMCIFIO, MMCIFParser, NeighborSearch, PDBIO, PDBParser, Superimposer
        from Bio.PDB.Polypeptide import is_aa
        from Bio.PDB.vectors import calc_dihedral

        bio_aligner = Align.PairwiseAligner()
        bio_aligner.mode = "global"
        bio_aligner.match_score = 2.0
        bio_aligner.mismatch_score = -1.0
        bio_aligner.open_gap_score = -1.0
        bio_aligner.extend_gap_score = -0.5

        def bio_read(path: str):
            parser = MMCIFParser(QUIET=True) if _is_cif(path) else PDBParser(QUIET=True)
            return parser.get_structure("model", path)

        def bio_write(obj: Any, out_path: str):
            io = MMCIFIO() if _is_cif(out_path) else PDBIO()
            io.set_structure(obj)
            io.save(out_path)

        def bio_select_within10_chainA(obj: Any) -> int:
            atoms = list(obj.get_atoms())
            chain_a_atoms = [
                atom for atom in atoms if atom.get_parent().get_parent().id == "A"
            ]
            if not chain_a_atoms:
                return 0
            search = NeighborSearch(atoms)
            selected = set()
            for atom in chain_a_atoms:
                for neighbor in search.search(atom.coord, 10.0, level="A"):
                    selected.add(id(neighbor))
            return len(selected)

        def bio_get_aa_seq_len(obj: Any) -> int:
            total = 0
            for model in obj:
                for chain in model:
                    for residue in chain:
                        if is_aa(residue, standard=False):
                            total += 1
            return total

        def bio_chain_seq_ca(obj: Any, chain_id: str) -> tuple[str, list[Any]]:
            seq = []
            atoms = []
            for atom in obj.get_atoms():
                if atom.get_name() != "CA":
                    continue
                residue = atom.get_parent()
                chain = residue.get_parent()
                if chain.id != chain_id:
                    continue
                seq.append(AA3_TO_1.get(residue.get_resname(), "X"))
                atoms.append(atom)
            return "".join(seq), atoms

        def bio_rmsd_ca_shift(reference_obj: Any, mobile_obj: Any) -> float:
            fixed_ca = [atom for atom in reference_obj.get_atoms() if atom.get_name() == "CA"]
            mobile_ca = [atom for atom in mobile_obj.get_atoms() if atom.get_name() == "CA"]
            n = min(len(fixed_ca), len(mobile_ca))
            if n == 0:
                return 0.0
            squared_sum = 0.0
            for fixed_atom, mobile_atom in zip(fixed_ca[:n], mobile_ca[:n]):
                distance = fixed_atom - mobile_atom
                squared_sum += distance * distance
            return (squared_sum / n) ** 0.5

        def bio_dihedral_ca(obj: Any) -> float:
            ca_atoms = [atom for atom in obj.get_atoms() if atom.get_name() == "CA"]
            if len(ca_atoms) < 4:
                return 0.0
            angle = calc_dihedral(
                ca_atoms[0].get_vector(),
                ca_atoms[1].get_vector(),
                ca_atoms[2].get_vector(),
                ca_atoms[3].get_vector(),
            )
            return float(np.degrees(angle))

        def bio_align_seq_chainA(reference_path: str, mobile_path: str) -> float:
            fixed = bio_read(reference_path)
            mobile = bio_read(mobile_path)
            fixed_seq, fixed_atoms = bio_chain_seq_ca(fixed, "A")
            mobile_seq, mobile_atoms = bio_chain_seq_ca(mobile, "A")
            if not fixed_atoms or not mobile_atoms:
                return 0.0
            alignment = bio_aligner.align(fixed_seq, mobile_seq)[0]
            fixed_sel = []
            mobile_sel = []
            for (fixed_start, fixed_end), (mobile_start, mobile_end) in zip(
                alignment.aligned[0], alignment.aligned[1]
            ):
                block_len = min(fixed_end - fixed_start, mobile_end - mobile_start)
                for offset in range(block_len):
                    fixed_sel.append(fixed_atoms[fixed_start + offset])
                    mobile_sel.append(mobile_atoms[mobile_start + offset])
            if not fixed_sel:
                return 0.0
            sup = Superimposer()
            sup.set_atoms(fixed_sel, mobile_sel)
            return float(sup.rms)

        def bio_align_ca_pair(reference_path: str, mobile_path: str) -> float:
            fixed = bio_read(reference_path)
            mobile = bio_read(mobile_path)
            fixed_ca = [atom for atom in fixed.get_atoms() if atom.get_name() == "CA"]
            mobile_ca = [atom for atom in mobile.get_atoms() if atom.get_name() == "CA"]
            n = min(len(fixed_ca), len(mobile_ca))
            if n == 0:
                return 0.0
            sup = Superimposer()
            sup.set_atoms(fixed_ca[:n], mobile_ca[:n])
            return float(sup.rms)

        def bio_hbond_count(obj: Any) -> int:
            # BioPython has no built-in H-bond function; use NeighborSearch-based D-A distance proxy
            from Bio.PDB import NeighborSearch
            atoms = list(obj.get_atoms())
            ns = NeighborSearch(atoms)
            count = 0
            for donor in atoms:
                if donor.element not in ("N", "O"):
                    continue
                for nb in ns.search(donor.coord, 3.5, level="A"):
                    if nb.element == "O" and nb != donor:
                        res_d = donor.get_parent()
                        res_a = nb.get_parent()
                        if res_d != res_a:
                            count += 1
            return count

        backends.append(
            Backend(
                name="biopython",
                read=bio_read,
                write=bio_write,
                select_within10_chainA=bio_select_within10_chainA,
                get_aa_seq_len=bio_get_aa_seq_len,
                rmsd_ca_shift=bio_rmsd_ca_shift,
                dihedral_ca=bio_dihedral_ca,
                align_seq_chainA=bio_align_seq_chainA,
                align_ca_pair=bio_align_ca_pair,
                hbond_count=bio_hbond_count,
            )
        )
    except Exception as exc:
        _report_backend_skip("biopython", exc)
        pass

    try:
        import biotite.structure as struc
        import biotite.structure.io.pdb as pdb
        import biotite.structure.io.pdbx as pdbx
        import biotite.sequence as seq
        import biotite.sequence.align as seq_align

        bt_matrix = seq_align.SubstitutionMatrix.std_protein_matrix()

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
            cell_list = struc.CellList(obj, 10.0)
            masks = cell_list.get_atoms(ref.coord, 10.0, as_mask=True)
            return int(np.count_nonzero(np.any(masks, axis=0)))

        def bt_get_aa_seq_len(obj: Any) -> int:
            _, res_names = struc.get_residues(obj)
            return sum(1 for name in res_names if name in AA3_TO_1)

        def bt_chain_seq_ca(obj: Any, chain_id: str) -> tuple[str, Any]:
            mask = (obj.atom_name == "CA") & (obj.chain_id == chain_id)
            if not np.any(mask):
                return "", obj[mask]
            seq = "".join(AA3_TO_1.get(str(r), "X") for r in obj.res_name[mask])
            return seq, obj[mask]

        def bt_rmsd_ca_shift(reference_obj: Any, mobile_obj: Any) -> float:
            fixed_ca = reference_obj[reference_obj.atom_name == "CA"]
            mobile_ca = mobile_obj[mobile_obj.atom_name == "CA"]
            n = min(fixed_ca.array_length(), mobile_ca.array_length())
            if n == 0:
                return 0.0
            return float(struc.rmsd(fixed_ca[:n], mobile_ca[:n]))

        def bt_dihedral_ca(obj: Any) -> float:
            ca = obj[obj.atom_name == "CA"]
            if ca.array_length() < 4:
                return 0.0
            angle = struc.dihedral(ca.coord[0], ca.coord[1], ca.coord[2], ca.coord[3])
            return float(np.degrees(angle))

        def bt_align_seq_chainA(reference_path: str, mobile_path: str) -> float:
            fixed = bt_read(reference_path)
            mobile = bt_read(mobile_path)
            fixed_seq, fixed_ca = bt_chain_seq_ca(fixed, "A")
            mobile_seq, mobile_ca = bt_chain_seq_ca(mobile, "A")
            if fixed_ca.array_length() == 0 or mobile_ca.array_length() == 0:
                return 0.0
            alignment = seq_align.align_optimal(
                seq.ProteinSequence(fixed_seq),
                seq.ProteinSequence(mobile_seq),
                bt_matrix,
                max_number=1,
            )[0]
            trace = alignment.trace
            trace = trace[(trace[:, 0] >= 0) & (trace[:, 1] >= 0)]
            if trace.shape[0] == 0:
                return 0.0
            fixed_sel = fixed_ca[trace[:, 0]]
            mobile_sel = mobile_ca[trace[:, 1]]
            fitted, _ = struc.superimpose(fixed_sel, mobile_sel)
            return float(struc.rmsd(fixed_sel, fitted))

        def bt_align_ca_pair(reference_path: str, mobile_path: str) -> float:
            fixed = bt_read(reference_path)
            mobile = bt_read(mobile_path)
            fixed_ca = fixed[fixed.atom_name == "CA"]
            mobile_ca = mobile[mobile.atom_name == "CA"]
            n = min(fixed_ca.array_length(), mobile_ca.array_length())
            if n == 0:
                return 0.0
            fitted, _ = struc.superimpose(fixed_ca[:n], mobile_ca[:n])
            return float(struc.rmsd(fixed_ca[:n], fitted))

        def bt_hbond_count(obj: Any) -> int:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                triplets = struc.hbond(obj, cutoff_dist=2.5, cutoff_angle=120)
            return int(len(triplets))

        backends.append(
            Backend(
                name="biotite",
                read=bt_read,
                write=bt_write,
                select_within10_chainA=bt_select_within10_chainA,
                get_aa_seq_len=bt_get_aa_seq_len,
                rmsd_ca_shift=bt_rmsd_ca_shift,
                dihedral_ca=bt_dihedral_ca,
                align_seq_chainA=bt_align_seq_chainA,
                align_ca_pair=bt_align_ca_pair,
                hbond_count=bt_hbond_count,
            )
        )
    except Exception as exc:
        _report_backend_skip("biotite", exc)
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
        "hbond",
    ]
    rows: list[dict[str, Any]] = []

    print("library\tfile\toperation\tmean_s\tstd_s\tspeedup_vs_pdb_cpp\tpdb_cpp_superior")

    for file_path in files:
        per_file_op_times: dict[tuple[str, str], float] = {}

        with tempfile.TemporaryDirectory(prefix="pdb_cpp_benchmark_") as temp_dir:
            transformed_path = Path(temp_dir) / f"{file_path.stem}_transformed{file_path.suffix}"
            _create_transformed_structure(file_path, transformed_path)
            for backend in backends:
                out_path = str(Path(temp_dir) / f"{backend.name}_{file_path.name}")
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
                    transformed_obj = backend.read(str(transformed_path))
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
                        lambda b=backend, reference=obj, mobile=transformed_obj: b.rmsd_ca_shift(reference, mobile),
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
                        lambda b=backend, current=obj: b.dihedral_ca(current),
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
                        lambda b=backend, reference=str(file_path), mobile=str(transformed_path): b.align_seq_chainA(reference, mobile),
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
                        lambda b=backend, reference=str(file_path), mobile=str(transformed_path): b.align_ca_pair(reference, mobile),
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

                    hbond_mean, hbond_std = benchmark_call(
                        lambda b=backend, current=obj: b.hbond_count(current),
                        args.warmup,
                        args.runs,
                    )
                    per_file_op_times[(backend.name, "hbond")] = hbond_mean
                    backend_key_list.append((backend.name, "hbond"))
                    rows.append(
                        {
                            "library": backend.name,
                            "file": str(file_path),
                            "operation": "hbond",
                            "mean_s": hbond_mean,
                            "std_s": hbond_std,
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
