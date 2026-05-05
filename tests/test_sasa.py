#!/usr/bin/env python3
# coding: utf-8

import math

import numpy as np
import pytest

from pdb_cpp import Coor, analysis
from pdb_cpp.analysis import sasa as analysis_sasa

from .datafiles import MMCIF_1Y0M, PDB_SASA_DIMER, PDB_SASA_SINGLE
from .sasa_reference_values import BURIED_SURFACE_REFERENCES, SASA_REFERENCES


DIMER_REFERENCE_ABS = 1.0
REAL_STRUCTURE_REFERENCE_ABS = 20.0


def _assert_area_breakdown(actual, expected, abs_tol):
    assert actual["total"] == pytest.approx(expected["total"], abs=abs_tol)
    assert actual["polar"] == pytest.approx(expected["polar"], abs=abs_tol)
    assert actual["apolar"] == pytest.approx(expected["apolar"], abs=abs_tol)


def _assert_residue_areas_match(actual_entries, expected_entries, abs_tol):
    assert len(actual_entries) == len(expected_entries)
    for actual, expected in zip(actual_entries, expected_entries):
        assert actual["chain"] == expected["chain"]
        assert actual["resid"] == expected["resid"]
        assert actual["resname"] == expected["resname"]
        assert actual["area"] == pytest.approx(expected["area"], abs=abs_tol)
        assert actual["polar_area"] == pytest.approx(expected["polar_area"], abs=abs_tol)
        assert actual["apolar_area"] == pytest.approx(expected["apolar_area"], abs=abs_tol)


def _assert_buried_residue_areas_match(actual_entries, expected_entries, abs_tol):
    assert len(actual_entries) == len(expected_entries)
    for actual, expected in zip(actual_entries, expected_entries):
        assert actual["partner"] == expected["partner"]
        assert actual["chain"] == expected["chain"]
        assert actual["resid"] == expected["resid"]
        assert actual["resname"] == expected["resname"]
        assert actual["isolated_area"] == pytest.approx(expected["isolated_area"], abs=abs_tol)
        assert actual["complex_area"] == pytest.approx(expected["complex_area"], abs=abs_tol)
        assert actual["buried_area"] == pytest.approx(expected["buried_area"], abs=abs_tol)
        assert actual["isolated_polar_area"] == pytest.approx(
            expected["isolated_polar_area"], abs=abs_tol
        )
        assert actual["complex_polar_area"] == pytest.approx(expected["complex_polar_area"], abs=abs_tol)
        assert actual["buried_polar_area"] == pytest.approx(expected["buried_polar_area"], abs=abs_tol)
        assert actual["isolated_apolar_area"] == pytest.approx(
            expected["isolated_apolar_area"], abs=abs_tol
        )
        assert actual["complex_apolar_area"] == pytest.approx(
            expected["complex_apolar_area"], abs=abs_tol
        )
        assert actual["buried_apolar_area"] == pytest.approx(expected["buried_apolar_area"], abs=abs_tol)


def test_single_atom_sasa_matches_analytic_area():
    coor = Coor(PDB_SASA_SINGLE)

    results = analysis_sasa.sasa(coor, n_points=960, by_atom=True)
    result = results[0]

    expected = 4.0 * math.pi * (1.7 + 1.4) ** 2
    assert len(results) == 1
    assert result["total"] == pytest.approx(expected, abs=1e-5)
    assert result["polar"] == pytest.approx(0.0, abs=1e-8)
    assert result["apolar"] == pytest.approx(expected, abs=1e-5)
    assert np.sum(result["atom_areas"]) == pytest.approx(result["total"], abs=1e-6)


def test_dimer_buried_surface_is_positive():
    coor = Coor(PDB_SASA_DIMER)
    expected = BURIED_SURFACE_REFERENCES["sasa_dimer_A_B"]

    metrics = analysis_sasa.buried_surface_area(coor, "chain A", "chain B", n_points=960)[0]

    assert metrics["receptor_sasa"] == pytest.approx(expected["receptor_sasa"], abs=DIMER_REFERENCE_ABS)
    assert metrics["receptor_polar_sasa"] == pytest.approx(
        expected["receptor_polar_sasa"], abs=DIMER_REFERENCE_ABS
    )
    assert metrics["receptor_apolar_sasa"] == pytest.approx(
        expected["receptor_apolar_sasa"], abs=DIMER_REFERENCE_ABS
    )
    assert metrics["ligand_sasa"] == pytest.approx(expected["ligand_sasa"], abs=DIMER_REFERENCE_ABS)
    assert metrics["ligand_polar_sasa"] == pytest.approx(
        expected["ligand_polar_sasa"], abs=DIMER_REFERENCE_ABS
    )
    assert metrics["ligand_apolar_sasa"] == pytest.approx(
        expected["ligand_apolar_sasa"], abs=DIMER_REFERENCE_ABS
    )
    assert metrics["complex_sasa"] == pytest.approx(expected["complex_sasa"], abs=DIMER_REFERENCE_ABS)
    assert metrics["complex_polar_sasa"] == pytest.approx(
        expected["complex_polar_sasa"], abs=DIMER_REFERENCE_ABS
    )
    assert metrics["complex_apolar_sasa"] == pytest.approx(
        expected["complex_apolar_sasa"], abs=DIMER_REFERENCE_ABS
    )
    assert metrics["buried_surface"] == pytest.approx(expected["buried_surface"], abs=DIMER_REFERENCE_ABS)
    assert metrics["buried_polar_surface"] == pytest.approx(
        expected["buried_polar_surface"], abs=DIMER_REFERENCE_ABS
    )
    assert metrics["buried_apolar_surface"] == pytest.approx(
        expected["buried_apolar_surface"], abs=DIMER_REFERENCE_ABS
    )
    assert metrics["interface_area"] == pytest.approx(expected["interface_area"], abs=DIMER_REFERENCE_ABS)
    assert metrics["interface_polar_area"] == pytest.approx(
        expected["interface_polar_area"], abs=DIMER_REFERENCE_ABS
    )
    assert metrics["interface_apolar_area"] == pytest.approx(
        expected["interface_apolar_area"], abs=DIMER_REFERENCE_ABS
    )


def test_atom_areas_sum_to_total_on_real_structure():
    coor = Coor(PDB_SASA_DIMER)
    expected = SASA_REFERENCES["sasa_dimer"]

    results = analysis_sasa.sasa(coor, n_points=480, by_atom=True)
    result = results[0]

    _assert_area_breakdown(result, expected, DIMER_REFERENCE_ABS)
    assert result["atom_areas"].shape == (coor.models[0].len,)
    assert np.sum(result["atom_areas"]) == pytest.approx(result["total"], abs=1e-5)


def test_model_sasa_can_return_residue_breakdown():
    coor = Coor(PDB_SASA_DIMER)
    expected = SASA_REFERENCES["sasa_dimer"]

    results = analysis_sasa.sasa(coor, n_points=480, by_residue=True)
    result = results[0]

    assert "atom_areas" not in result
    _assert_area_breakdown(result, expected, DIMER_REFERENCE_ABS)
    _assert_residue_areas_match(result["residue_areas"], expected["residue_areas"], DIMER_REFERENCE_ABS)


def test_buried_surface_can_return_residue_decomposition():
    coor = Coor(PDB_SASA_DIMER)
    expected = BURIED_SURFACE_REFERENCES["sasa_dimer_A_B"]

    metrics = analysis_sasa.buried_surface_area(
        coor, "chain A", "chain B", n_points=480, by_residue=True
    )[0]

    _assert_residue_areas_match(
        metrics["receptor_residue_sasa"], expected["receptor_residue_sasa"], DIMER_REFERENCE_ABS
    )
    _assert_residue_areas_match(
        metrics["ligand_residue_sasa"], expected["ligand_residue_sasa"], DIMER_REFERENCE_ABS
    )
    _assert_residue_areas_match(
        metrics["complex_residue_sasa"], expected["complex_residue_sasa"], DIMER_REFERENCE_ABS
    )
    _assert_buried_residue_areas_match(
        metrics["residue_buried_surface"], expected["residue_buried_surface"], DIMER_REFERENCE_ABS
    )


def test_sasa_selection_returns_one_result_per_model():
    coor = Coor(PDB_SASA_DIMER)
    expected = SASA_REFERENCES["sasa_dimer_chain_A"]

    results = analysis_sasa.sasa(coor, selection="chain A", n_points=480, by_residue=True)

    assert len(results) == coor.model_num
    _assert_area_breakdown(results[0], expected, DIMER_REFERENCE_ABS)
    _assert_residue_areas_match(results[0]["residue_areas"], expected["residue_areas"], DIMER_REFERENCE_ABS)


def test_real_structure_sasa_reports_nonzero_polar_and_apolar_components():
    coor = Coor(MMCIF_1Y0M)
    expected = SASA_REFERENCES["1y0m"]

    result = analysis_sasa.sasa(coor, n_points=480)[0]

    _assert_area_breakdown(result, expected, REAL_STRUCTURE_REFERENCE_ABS)


def test_analysis_namespace_sasa_is_still_callable():
    coor = Coor(PDB_SASA_SINGLE)

    result = analysis.sasa(coor)[0]

    assert result["total"] == pytest.approx(result["polar"] + result["apolar"], abs=1e-5)
    assert result["apolar"] > 0.0


def test_analysis_sasa_module_imports_work():
    coor = Coor(PDB_SASA_SINGLE)

    result = analysis_sasa.sasa(coor)[0]

    assert result["total"] == pytest.approx(result["polar"] + result["apolar"], abs=1e-5)


def test_shape_complementarity_simple_dimer_is_positive_and_symmetric():
    coor = Coor(PDB_SASA_DIMER)

    forward = analysis_sasa.shape_complementarity(
        coor,
        "chain A",
        "chain B",
        dots_per_sq_angstrom=6.0,
        search_radius=1.0,
    )[0]
    reverse = analysis_sasa.shape_complementarity(
        coor,
        "chain B",
        "chain A",
        dots_per_sq_angstrom=6.0,
        search_radius=1.0,
    )[0]

    assert -1.0 <= forward["shape_complementarity"] <= 1.0
    assert forward["shape_complementarity"] == pytest.approx(
        reverse["shape_complementarity"],
        abs=1e-6,
    )
    assert forward["interface_dot_pairs"] > 0
    assert 0.0 <= forward["mean_interface_distance"] <= 1.0


def test_shape_complementarity_flat_analysis_alias_matches_module():
    coor = Coor(PDB_SASA_DIMER)

    flat_result = analysis.shape_complementarity(
        coor,
        "chain A",
        "chain B",
        dots_per_sq_angstrom=6.0,
        search_radius=1.0,
    )[0]
    module_result = analysis_sasa.shape_complementarity(
        coor,
        "chain A",
        "chain B",
        dots_per_sq_angstrom=6.0,
        search_radius=1.0,
    )[0]

    assert flat_result["shape_complementarity"] == pytest.approx(
        module_result["shape_complementarity"],
        abs=1e-6,
    )


def test_shape_complementarity_raises_without_interface_pairs():
    coor = Coor(PDB_SASA_DIMER)

    with pytest.raises(ValueError, match="No opposing interface surface dots"):
        analysis_sasa.shape_complementarity(
            coor,
            "chain A",
            "chain B",
            dots_per_sq_angstrom=6.0,
            search_radius=0.01,
        )