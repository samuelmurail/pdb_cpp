#!/usr/bin/env python3
# coding: utf-8

import math
import importlib.util

import numpy as np
import pytest

from pdb_cpp import Coor, analysis

from .datafiles import PDB_SASA_DIMER, PDB_SASA_SINGLE


def test_single_atom_sasa_matches_analytic_area():
    coor = Coor(PDB_SASA_SINGLE)

    results = analysis.sasa(coor, n_points=960, by_atom=True)
    result = results[0]

    expected = 4.0 * math.pi * (1.7 + 1.4) ** 2
    assert len(results) == 1
    assert result["total"] == pytest.approx(expected, abs=1e-5)
    assert np.sum(result["atom_areas"]) == pytest.approx(result["total"], abs=1e-6)


def test_dimer_buried_surface_is_positive():
    coor = Coor(PDB_SASA_DIMER)

    metrics = analysis.buried_surface_area(coor, "chain A", "chain B", n_points=960)[0]

    assert metrics["receptor_sasa"] > metrics["complex_sasa"] / 2.0
    assert metrics["ligand_sasa"] > metrics["complex_sasa"] / 2.0
    assert metrics["buried_surface"] > 0.0
    assert metrics["interface_area"] == pytest.approx(metrics["buried_surface"] / 2.0, abs=1e-6)


def test_atom_areas_sum_to_total_on_real_structure():
    coor = Coor(PDB_SASA_DIMER)

    results = analysis.sasa(coor, n_points=480, by_atom=True)
    result = results[0]

    assert result["atom_areas"].shape == (coor.models[0].len,)
    assert np.sum(result["atom_areas"]) == pytest.approx(result["total"], abs=1e-5)


def test_model_sasa_can_return_residue_breakdown():
    coor = Coor(PDB_SASA_DIMER)

    results = analysis.sasa(coor, n_points=480, by_residue=True)
    result = results[0]

    assert "atom_areas" not in result
    assert len(result["residue_areas"]) == 2
    assert sum(entry["area"] for entry in result["residue_areas"]) == pytest.approx(
        result["total"],
        abs=1e-5,
    )
    assert result["residue_areas"][0]["chain"] == "A"
    assert result["residue_areas"][1]["chain"] == "B"


def test_buried_surface_can_return_residue_decomposition():
    coor = Coor(PDB_SASA_DIMER)

    metrics = analysis.buried_surface_area(coor, "chain A", "chain B", n_points=480, by_residue=True)[0]

    assert len(metrics["receptor_residue_sasa"]) == 1
    assert len(metrics["ligand_residue_sasa"]) == 1
    assert len(metrics["complex_residue_sasa"]) == 2
    assert len(metrics["residue_buried_surface"]) == 2
    assert sum(entry["buried_area"] for entry in metrics["residue_buried_surface"]) == pytest.approx(
        metrics["buried_surface"],
        abs=1e-5,
    )
    assert all(entry["buried_area"] > 0.0 for entry in metrics["residue_buried_surface"])


def test_invalid_sasa_backend_raises():
    coor = Coor(PDB_SASA_DIMER)

    with pytest.raises(ValueError, match="backend"):
        analysis.buried_surface_area(coor, "chain A", "chain B", backend="bad-backend")


@pytest.mark.skipif(importlib.util.find_spec("freesasa") is None, reason="freesasa not installed")
def test_freesasa_backend_matches_direct_freesasa_totals(tmp_path):
    import freesasa

    coor = Coor(PDB_SASA_DIMER)
    metrics = analysis.buried_surface_area(coor, "chain A", "chain B", backend="freesasa", by_residue=True)[0]

    rec = tmp_path / "rec.pdb"
    lig = tmp_path / "lig.pdb"
    com = tmp_path / "com.pdb"
    coor.select_atoms("chain A").write(str(rec))
    coor.select_atoms("chain B").write(str(lig))
    coor.select_atoms("(chain A) or (chain B)").write(str(com))

    rec_total = freesasa.calc(freesasa.Structure(str(rec))).totalArea()
    lig_total = freesasa.calc(freesasa.Structure(str(lig))).totalArea()
    com_total = freesasa.calc(freesasa.Structure(str(com))).totalArea()

    assert metrics["receptor_sasa"] == pytest.approx(rec_total, abs=1e-8)
    assert metrics["ligand_sasa"] == pytest.approx(lig_total, abs=1e-8)
    assert metrics["complex_sasa"] == pytest.approx(com_total, abs=1e-8)
    assert sum(entry["buried_area"] for entry in metrics["residue_buried_surface"]) == pytest.approx(
        metrics["buried_surface"],
        abs=1e-8,
    )


def test_sasa_selection_returns_one_result_per_model():
    coor = Coor(PDB_SASA_DIMER)

    results = analysis.sasa(coor, selection="chain A", by_residue=True)

    assert len(results) == coor.model_num
    assert len(results[0]["residue_areas"]) == 1
    assert results[0]["residue_areas"][0]["chain"] == "A"