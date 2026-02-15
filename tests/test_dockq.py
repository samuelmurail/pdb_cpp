#!/usr/bin/env python3
# coding: utf-8

"""DockQ reference-value regression tests."""

from .datafiles import DOCKQ_MODEL, DOCKQ_NATIVE, PDB_1JD4, PDB_1RXZ, PDB_1RXZ_Colabfold, PDB_5M6N
from .dockq_reference_values import DOCKQ_REFERENCES
from pdb_cpp import Coor, analysis
import pytest


def test_dockq_reference_1rxz_colabfold_vs_native():
    model_coor = Coor(PDB_1RXZ_Colabfold)
    native_coor = Coor(PDB_1RXZ)

    result = analysis.dockQ(model_coor, native_coor)
    ref = DOCKQ_REFERENCES["1rxz_colabfold_vs_1rxz"]

    assert result["DockQ"][0] == pytest.approx(ref["DockQ"], abs=0.001)
    assert result["Fnat"][0] == pytest.approx(ref["Fnat"], abs=0.001)
    assert result["Fnonnat"][0] == pytest.approx(ref["Fnonnat"], abs=0.001)
    assert result["iRMS"][0] == pytest.approx(ref["iRMS"], abs=0.001)
    assert result["LRMS"][0] == pytest.approx(ref["LRMS"], abs=0.001)


def test_dockq_reference_model_vs_native():
    model_coor = Coor(DOCKQ_MODEL)
    native_coor = Coor(DOCKQ_NATIVE)

    result = analysis.dockQ(model_coor, native_coor)
    ref = DOCKQ_REFERENCES["model_vs_native"]

    assert result["DockQ"][0] == pytest.approx(ref["DockQ"], abs=0.001)
    assert result["Fnat"][0] == pytest.approx(ref["Fnat"], abs=0.001)
    assert result["Fnonnat"][0] == pytest.approx(ref["Fnonnat"], abs=0.001)
    assert result["iRMS"][0] == pytest.approx(ref["iRMS"], abs=0.001)
    assert result["LRMS"][0] == pytest.approx(ref["LRMS"], abs=0.001)


def test_dockq_reference_1jd4_vs_5m6n_score_target():
    model_coor = Coor(PDB_1JD4)
    native_coor = Coor(PDB_5M6N)

    result = analysis.dockQ(model_coor, native_coor)
    ref = DOCKQ_REFERENCES["1jd4_vs_5m6n"]

    assert result["DockQ"][0] == pytest.approx(ref["DockQ"], abs=0.005)
    assert result["Fnat"][0] == pytest.approx(ref["Fnat"], abs=0.001)
    assert result["Fnonnat"][0] == pytest.approx(ref["Fnonnat"], abs=0.001)