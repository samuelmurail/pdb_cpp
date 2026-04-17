#!/usr/bin/env python3
# coding: utf-8

import math

import numpy as np
import pytest

from pdb_cpp import Coor, analysis, sasa

from .datafiles import MMCIF_1Y0M, PDB_SASA_DIMER, PDB_SASA_SINGLE


def test_single_atom_sasa_matches_analytic_area():
    coor = Coor(PDB_SASA_SINGLE)

    results = sasa.sasa(coor, n_points=960, by_atom=True)
    result = results[0]

    expected = 4.0 * math.pi * (1.7 + 1.4) ** 2
    assert len(results) == 1
    assert result["total"] == pytest.approx(expected, abs=1e-5)
    assert result["polar"] == pytest.approx(0.0, abs=1e-8)
    assert result["apolar"] == pytest.approx(expected, abs=1e-5)
    assert np.sum(result["atom_areas"]) == pytest.approx(result["total"], abs=1e-6)


def test_dimer_buried_surface_is_positive():
    coor = Coor(PDB_SASA_DIMER)

    metrics = sasa.buried_surface_area(coor, "chain A", "chain B", n_points=960)[0]

    assert metrics["receptor_sasa"] > metrics["complex_sasa"] / 2.0
    assert metrics["ligand_sasa"] > metrics["complex_sasa"] / 2.0
    assert metrics["buried_surface"] > 0.0
    assert metrics["receptor_sasa"] == pytest.approx(
        metrics["receptor_polar_sasa"] + metrics["receptor_apolar_sasa"], abs=1e-6
    )
    assert metrics["ligand_sasa"] == pytest.approx(
        metrics["ligand_polar_sasa"] + metrics["ligand_apolar_sasa"], abs=1e-6
    )
    assert metrics["complex_sasa"] == pytest.approx(
        metrics["complex_polar_sasa"] + metrics["complex_apolar_sasa"], abs=1e-6
    )
    assert metrics["buried_surface"] == pytest.approx(
        metrics["buried_polar_surface"] + metrics["buried_apolar_surface"], abs=1e-6
    )
    assert metrics["interface_area"] == pytest.approx(metrics["buried_surface"] / 2.0, abs=1e-6)
    assert metrics["interface_area"] == pytest.approx(
        metrics["interface_polar_area"] + metrics["interface_apolar_area"], abs=1e-6
    )


def test_atom_areas_sum_to_total_on_real_structure():
    coor = Coor(PDB_SASA_DIMER)

    results = sasa.sasa(coor, n_points=480, by_atom=True)
    result = results[0]

    assert result["atom_areas"].shape == (coor.models[0].len,)
    assert np.sum(result["atom_areas"]) == pytest.approx(result["total"], abs=1e-5)


def test_model_sasa_can_return_residue_breakdown():
    coor = Coor(PDB_SASA_DIMER)

    results = sasa.sasa(coor, n_points=480, by_residue=True)
    result = results[0]

    assert "atom_areas" not in result
    assert len(result["residue_areas"]) == 2
    assert sum(entry["area"] for entry in result["residue_areas"]) == pytest.approx(
        result["total"],
        abs=1e-5,
    )
    assert all(
        entry["area"] == pytest.approx(entry["polar_area"] + entry["apolar_area"], abs=1e-6)
        for entry in result["residue_areas"]
    )
    assert result["residue_areas"][0]["chain"] == "A"
    assert result["residue_areas"][1]["chain"] == "B"


def test_buried_surface_can_return_residue_decomposition():
    coor = Coor(PDB_SASA_DIMER)

    metrics = sasa.buried_surface_area(coor, "chain A", "chain B", n_points=480, by_residue=True)[0]

    assert len(metrics["receptor_residue_sasa"]) == 1
    assert len(metrics["ligand_residue_sasa"]) == 1
    assert len(metrics["complex_residue_sasa"]) == 2
    assert len(metrics["residue_buried_surface"]) == 2
    assert sum(entry["buried_area"] for entry in metrics["residue_buried_surface"]) == pytest.approx(
        metrics["buried_surface"],
        abs=1e-5,
    )
    assert sum(entry["buried_polar_area"] for entry in metrics["residue_buried_surface"]) == pytest.approx(
        metrics["buried_polar_surface"],
        abs=1e-5,
    )
    assert sum(entry["buried_apolar_area"] for entry in metrics["residue_buried_surface"]) == pytest.approx(
        metrics["buried_apolar_surface"],
        abs=1e-5,
    )
    assert all(entry["buried_area"] > 0.0 for entry in metrics["residue_buried_surface"])


def test_sasa_selection_returns_one_result_per_model():
    coor = Coor(PDB_SASA_DIMER)

    results = sasa.sasa(coor, selection="chain A", by_residue=True)

    assert len(results) == coor.model_num
    assert len(results[0]["residue_areas"]) == 1
    assert results[0]["residue_areas"][0]["chain"] == "A"


def test_real_structure_sasa_reports_nonzero_polar_and_apolar_components():
    coor = Coor(MMCIF_1Y0M)

    result = sasa.sasa(coor, n_points=480)[0]

    assert result["total"] == pytest.approx(result["polar"] + result["apolar"], abs=1e-5)
    assert result["polar"] > 0.0
    assert result["apolar"] > 0.0


def test_analysis_sasa_backwards_compatibility():
    coor = Coor(PDB_SASA_SINGLE)

    result = analysis.sasa(coor)[0]

    assert result["total"] == pytest.approx(result["polar"] + result["apolar"], abs=1e-5)
    assert result["apolar"] > 0.0