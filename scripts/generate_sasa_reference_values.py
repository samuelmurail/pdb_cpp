#!/usr/bin/env python3
# coding: utf-8

"""Generate SASA regression references from the FreeSASA CLI.

The reference values are produced by:
1. rewriting PDB occupancy values to the same element radii used by
   ``pdb_cpp.analysis.sasa``
2. running ``freesasa --radius-from-occupancy --hetatm --format=json``
3. writing the resulting values to ``tests/sasa_reference_values.py``

This script expects the ``freesasa`` CLI to be available, either directly on
``PATH`` or through ``conda run -n <env> freesasa``.
"""

from __future__ import annotations

import argparse
import ast
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from pprint import pformat


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = PROJECT_ROOT / "tests" / "sasa_reference_values.py"
TEST_INPUTS = PROJECT_ROOT / "tests" / "input"

_VDW_RADII = {
    "H": 1.10,
    "D": 1.10,
    "HE": 1.40,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "NE": 1.54,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "AR": 1.88,
    "SE": 1.90,
    "BR": 1.85,
    "KR": 2.02,
    "I": 1.98,
    "MG": 1.73,
    "NA": 2.27,
    "K": 2.75,
    "CA": 2.31,
    "MN": 1.73,
    "FE": 1.72,
    "CU": 1.40,
    "ZN": 1.39,
}

_TWO_LETTER_ELEMENTS = {"CL", "BR", "SE", "NA", "MG", "ZN", "FE", "CA", "MN", "CU"}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Python file to write. Default: tests/sasa_reference_values.py",
    )
    parser.add_argument(
        "--freesasa-env",
        default="freesasa",
        help="Conda environment that provides the freesasa CLI. Use '' to call freesasa directly.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Validate that the generated content matches the existing output file.",
    )
    return parser.parse_args()


def _infer_element(atom_name: str, elem: str) -> str:
    symbol = elem.strip().upper()
    if symbol:
        return symbol
    letters = "".join(char for char in atom_name.upper() if char.isalpha())
    if not letters:
        return "C"
    if len(letters) >= 2 and letters[:2] in _TWO_LETTER_ELEMENTS:
        return letters[:2]
    return letters[0]


def _rewrite_pdb_with_radii(source_path: Path, destination_path: Path, chain: str | None = None) -> None:
    with source_path.open() as handle, destination_path.open("w") as output_handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                if chain is not None and line[21].strip() != chain:
                    continue
                atom_name = line[12:16]
                element = line[76:78]
                radius = _VDW_RADII.get(_infer_element(atom_name, element), 1.70)
                line = f"{line[:54]}{radius:6.2f}{line[60:]}"
            output_handle.write(line)


def _run_freesasa_json(structure_path: Path, env_name: str) -> dict:
    if env_name:
        command = ["conda", "run", "-n", env_name, "freesasa"]
    else:
        command = ["freesasa"]
    command.extend(
        [
            "--format=json",
            "--depth=residue",
            "--radius-from-occupancy",
            "--hetatm",
            str(structure_path),
        ]
    )
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    return json.loads(completed.stdout)


def _structure_reference(source_path: Path, env_name: str, chain: str | None = None) -> dict:
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as handle:
        temp_path = Path(handle.name)
    try:
        _rewrite_pdb_with_radii(source_path, temp_path, chain=chain)
        result = _run_freesasa_json(temp_path, env_name)["results"][0]["structure"][0]
    finally:
        temp_path.unlink(missing_ok=True)

    residue_areas = []
    for chain_entry in result["chains"]:
        for residue in chain_entry["residues"]:
            residue_areas.append(
                {
                    "chain": chain_entry["label"],
                    "resid": int(residue["number"]),
                    "resname": residue["name"],
                    "area": residue["area"]["total"],
                    "polar_area": residue["area"]["polar"],
                    "apolar_area": residue["area"]["apolar"],
                }
            )

    structure_area = result["area"]
    reference = {
        "total": structure_area["total"],
        "polar": structure_area["polar"],
        "apolar": structure_area["apolar"],
    }
    if residue_areas:
        reference["residue_areas"] = residue_areas
    return reference


def _buried_surface_reference(complex_reference: dict, receptor_reference: dict, ligand_reference: dict) -> dict:
    reference = {
        "receptor_sasa": receptor_reference["total"],
        "receptor_polar_sasa": receptor_reference["polar"],
        "receptor_apolar_sasa": receptor_reference["apolar"],
        "ligand_sasa": ligand_reference["total"],
        "ligand_polar_sasa": ligand_reference["polar"],
        "ligand_apolar_sasa": ligand_reference["apolar"],
        "complex_sasa": complex_reference["total"],
        "complex_polar_sasa": complex_reference["polar"],
        "complex_apolar_sasa": complex_reference["apolar"],
    }
    reference["buried_surface"] = (
        reference["receptor_sasa"] + reference["ligand_sasa"] - reference["complex_sasa"]
    )
    reference["buried_polar_surface"] = (
        reference["receptor_polar_sasa"]
        + reference["ligand_polar_sasa"]
        - reference["complex_polar_sasa"]
    )
    reference["buried_apolar_surface"] = (
        reference["receptor_apolar_sasa"]
        + reference["ligand_apolar_sasa"]
        - reference["complex_apolar_sasa"]
    )
    reference["interface_area"] = reference["buried_surface"] / 2.0
    reference["interface_polar_area"] = reference["buried_polar_surface"] / 2.0
    reference["interface_apolar_area"] = reference["buried_apolar_surface"] / 2.0
    reference["receptor_residue_sasa"] = receptor_reference["residue_areas"]
    reference["ligand_residue_sasa"] = ligand_reference["residue_areas"]
    reference["complex_residue_sasa"] = complex_reference["residue_areas"]

    complex_lookup = {
        (entry["chain"], entry["resid"], entry["resname"]): entry
        for entry in complex_reference["residue_areas"]
    }
    residue_buried_surface = []
    for partner_name, isolated_entries in (
        ("receptor", receptor_reference["residue_areas"]),
        ("ligand", ligand_reference["residue_areas"]),
    ):
        for entry in isolated_entries:
            complex_entry = complex_lookup[(entry["chain"], entry["resid"], entry["resname"])]
            residue_buried_surface.append(
                {
                    "partner": partner_name,
                    "chain": entry["chain"],
                    "resid": entry["resid"],
                    "resname": entry["resname"],
                    "isolated_area": entry["area"],
                    "complex_area": complex_entry["area"],
                    "buried_area": entry["area"] - complex_entry["area"],
                    "isolated_polar_area": entry["polar_area"],
                    "complex_polar_area": complex_entry["polar_area"],
                    "buried_polar_area": entry["polar_area"] - complex_entry["polar_area"],
                    "isolated_apolar_area": entry["apolar_area"],
                    "complex_apolar_area": complex_entry["apolar_area"],
                    "buried_apolar_area": entry["apolar_area"] - complex_entry["apolar_area"],
                }
            )
    reference["residue_buried_surface"] = residue_buried_surface
    return reference


def _generate_reference_payload(env_name: str) -> tuple[dict, dict]:
    dimer_path = TEST_INPUTS / "sasa_dimer.pdb"
    y0m_path = TEST_INPUTS / "1y0m.pdb"

    dimer_complex = _structure_reference(dimer_path, env_name)
    dimer_chain_a = _structure_reference(dimer_path, env_name, chain="A")
    dimer_chain_b = _structure_reference(dimer_path, env_name, chain="B")
    y0m_reference = _structure_reference(y0m_path, env_name)

    sasa_references = {
        "1y0m": {
            "total": y0m_reference["total"],
            "polar": y0m_reference["polar"],
            "apolar": y0m_reference["apolar"],
        },
        "sasa_dimer": dimer_complex,
        "sasa_dimer_chain_A": dimer_chain_a,
    }
    buried_surface_references = {
        "sasa_dimer_A_B": _buried_surface_reference(dimer_complex, dimer_chain_a, dimer_chain_b)
    }
    return sasa_references, buried_surface_references


def _build_output_text(sasa_references: dict, buried_surface_references: dict) -> str:
    header = '''#!/usr/bin/env python3
# coding: utf-8

"""Reference SASA outputs captured from FreeSASA CLI.

These values were generated by ``scripts/generate_sasa_reference_values.py``
with FreeSASA using ``--radius-from-occupancy --hetatm`` after rewriting the
PDB occupancy column to the element-based radii defined in
``pdb_cpp.analysis.sasa``.
"""

'''
    body = (
        f"SASA_REFERENCES = {pformat(sasa_references, sort_dicts=False)}\n\n\n"
        f"BURIED_SURFACE_REFERENCES = {pformat(buried_surface_references, sort_dicts=False)}\n"
    )
    return header + body


def _assert_ast_equal(existing_text: str, generated_text: str) -> None:
    existing_tree = ast.dump(ast.parse(existing_text), include_attributes=False)
    generated_tree = ast.dump(ast.parse(generated_text), include_attributes=False)
    if existing_tree != generated_tree:
        raise SystemExit("Generated SASA references differ from the checked-in file.")


def main() -> int:
    args = _parse_args()
    sasa_references, buried_surface_references = _generate_reference_payload(args.freesasa_env)
    output_text = _build_output_text(sasa_references, buried_surface_references)

    if args.check:
        existing_text = args.output.read_text()
        _assert_ast_equal(existing_text, output_text)
        print(f"SASA references are up to date: {args.output}")
        return 0

    args.output.write_text(output_text)
    print(f"Wrote SASA references to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())