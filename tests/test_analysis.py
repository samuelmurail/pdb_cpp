#!/usr/bin/env python3
# coding: utf-8

"""RMSD and DockQ tests."""

from .datafiles import (
    DOCKQ_MODEL,
    DOCKQ_NATIVE,
    PDB_1JD4,
    PDB_1RXZ,
    PDB_1RXZ_Colabfold,
    PDB_1U85,
    PDB_1UBD,
    PDB_5M6N,
)
from pdb_cpp import Coor, alignment, analysis, core
import pytest


def test_measure_rmsd(capsys):
    coor_1 = Coor(PDB_1U85)
    coor_2 = Coor(PDB_1UBD)

    seq_1 = coor_1.get_aa_seq()
    seq_2 = coor_2.get_aa_seq()
    align_seq_1, align_seq_2, _ = alignment.align_seq(seq_1["A"], seq_2["C"])
    alignment.print_align_seq(align_seq_1, align_seq_2, line_len=80)
    _ = capsys.readouterr()

    index_1, index_2 = core.get_common_atoms(
        coor_1, coor_2, chain_1=["A"], chain_2=["C"]
    )

    assert len(index_1) == len(index_2) == 112
    rmsd = analysis.rmsd(coor_1, coor_2, index_list=[index_1, index_2])
    assert rmsd[0] == pytest.approx(70.38518415577853, 0.0001)

    core.coor_align(coor_1, coor_2, index_1, index_2, frame_ref=0)

    rmsds = analysis.rmsd(coor_1, coor_2, index_list=[index_1, index_2])
    expected_rmsds = [
        5.1201007697145995,
        4.325464568500979,
        3.814838140492011,
        3.7162291711703648,
        3.885813512555148,
        5.148095052210754,
        5.296391465950272,
        4.135615244634669,
        3.8189144358192806,
        4.597449831608669,
        5.271310413581032,
        5.517576912040033,
        4.6082437633178115,
        4.2097575131149885,
        4.996842582024358,
        5.006402154252272,
        5.256112097498127,
        3.7419617535551613,
        4.184792438296149,
        4.178818177627158,
    ]

    for expected_rmsd, rmsd_val in zip(expected_rmsds, rmsds):
        assert expected_rmsd == pytest.approx(rmsd_val, 0.0001)


def test_dockq_bad():
    model_coor = Coor(PDB_1JD4)
    native_coor = Coor(PDB_5M6N)

    dockq = analysis.dockQ(model_coor, native_coor)
    assert dockq["Fnat"][0] == 0.0
    assert dockq["Fnonnat"][0] == 1.0
    assert pytest.approx(dockq["LRMS"][0], 0.1) == 54.0
    assert pytest.approx(dockq["iRMS"][0], 0.5) == 15.631
    assert pytest.approx(dockq["DockQ"][0], 0.5) == 0.010


def test_dockq_good():
    model_coor = Coor(PDB_1RXZ_Colabfold)
    native_coor = Coor(PDB_1RXZ)

    dockq = analysis.dockQ(model_coor, native_coor)

    assert pytest.approx(dockq["DockQ"][0], 0.01) == 0.934
    assert pytest.approx(dockq["Fnat"][0], 0.01) == 0.963
    assert pytest.approx(dockq["Fnonnat"][0], 0.01) == 0.088
    assert pytest.approx(dockq["LRMS"][0], 0.1) == 1.050
    assert pytest.approx(dockq["iRMS"][0], 0.5) == 0.618


def test_dockq_model():
    model_coor = Coor(DOCKQ_MODEL)
    native_coor = Coor(DOCKQ_NATIVE)

    dockq = analysis.dockQ(model_coor, native_coor)

    assert pytest.approx(dockq["DockQ"][0], 0.5) == 0.7
    assert pytest.approx(dockq["Fnat"][0], 0.01) == 0.533
    assert pytest.approx(dockq["Fnonnat"][0], 0.01) == 0.238
    assert pytest.approx(dockq["LRMS"][0], 0.1) == 1.516
    assert pytest.approx(dockq["iRMS"][0], 0.5) == 1.232
