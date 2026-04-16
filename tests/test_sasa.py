#!/usr/bin/env python3
# coding: utf-8

import math

import numpy as np
import pytest

from pdb_cpp import Coor, sasa

from .datafiles import PDB_SASA_DIMER, PDB_SASA_SINGLE


def test_single_atom_sasa_matches_analytic_area():
    coor = Coor(PDB_SASA_SINGLE)

    result = coor.models[0].sasa(n_points=960, by_atom=True)

    expected = 4.0 * math.pi * (1.7 + 1.4) ** 2
    assert result["total"] == pytest.approx(expected, abs=1e-5)
    assert np.sum(result["atom_areas"]) == pytest.approx(result["total"], abs=1e-6)


def test_dimer_buried_surface_is_positive():
    coor = Coor(PDB_SASA_DIMER)

    metrics = sasa.buried_surface_area(coor, "chain A", "chain B", n_points=960)

    assert metrics["receptor_sasa"] > metrics["complex_sasa"] / 2.0
    assert metrics["ligand_sasa"] > metrics["complex_sasa"] / 2.0
    assert metrics["buried_surface"] > 0.0
    assert metrics["interface_area"] == pytest.approx(metrics["buried_surface"] / 2.0, abs=1e-6)


def test_atom_areas_sum_to_total_on_real_structure():
    coor = Coor(PDB_SASA_DIMER)

    result = coor.models[0].sasa(n_points=480, by_atom=True)

    assert result["atom_areas"].shape == (coor.models[0].len,)
    assert np.sum(result["atom_areas"]) == pytest.approx(result["total"], abs=1e-5)


def test_model_sasa_can_return_residue_breakdown():
    coor = Coor(PDB_SASA_DIMER)

    result = coor.models[0].sasa(n_points=480, by_residue=True)

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

    metrics = sasa.buried_surface_area(coor, "chain A", "chain B", n_points=480, by_residue=True)

    assert len(metrics["receptor_residue_sasa"]) == 1
    assert len(metrics["ligand_residue_sasa"]) == 1
    assert len(metrics["complex_residue_sasa"]) == 2
    assert len(metrics["residue_buried_surface"]) == 2
    assert sum(entry["buried_area"] for entry in metrics["residue_buried_surface"]) == pytest.approx(
        metrics["buried_surface"],
        abs=1e-5,
    )
    assert all(entry["buried_area"] > 0.0 for entry in metrics["residue_buried_surface"])