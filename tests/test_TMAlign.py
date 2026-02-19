#!/usr/bin/env python3
# coding: utf-8

"""TM-align and secondary-structure regression tests."""

from __future__ import annotations

import pytest

from .datafiles import (
    DOCKQ_MODEL,
    DOCKQ_NATIVE,
    PDB_1JD4,
    PDB_1RXZ,
    PDB_1RXZ_Colabfold,
    PDB_1UBD,
    PDB_1U85,
    PDB_1Y0M,
    PDB_2RRI,
    PDB_5M6N,
)
from pdb_cpp import Coor, TMalign
from pdb_cpp.core import tmalign_ca


def test_tmalign_matches_usalign_self_alignment():
    coor = Coor(PDB_1Y0M)

    result = tmalign_ca(coor, coor, chain_1=["A"], chain_2=["A"], mm=1)

    assert result.L_ali == 61
    assert result.rmsd == pytest.approx(0.0, abs=1e-6)
    assert sorted([result.TM1, result.TM2]) == pytest.approx([1.0, 1.0], abs=1e-6)


def test_tmalign_matches_usalign_nontrivial_pair():
    coor_1 = Coor(PDB_1Y0M)
    coor_2 = Coor(PDB_1UBD)

    result = tmalign_ca(coor_1, coor_2, chain_1=["A"], chain_2=["C"], mm=1)

    assert result.L_ali == 26
    assert result.rmsd == pytest.approx(3.62, abs=5e-2)
    assert sorted([result.TM1, result.TM2]) == pytest.approx(
        sorted([0.21408, 0.14362]), abs=2e-4
    )


def test_tmalign_matches_usalign_multimodel_first_frame():
    coor_1 = Coor(PDB_1Y0M)
    coor_2 = Coor(PDB_2RRI)

    result = tmalign_ca(coor_1, coor_2, chain_1=["A"], chain_2=["A"], mm=1)

    assert result.L_ali == 11
    assert result.rmsd == pytest.approx(1.88, abs=5e-2)
    assert sorted([result.TM1, result.TM2]) == pytest.approx(
        sorted([0.17387, 0.13299]), abs=2e-4
    )


def test_tmalign_alignment_strings_consistent():
    coor_1 = Coor(PDB_1Y0M)
    coor_2 = Coor(PDB_2RRI)

    result = tmalign_ca(coor_1, coor_2, chain_1=["A"], chain_2=["A"], mm=1)

    assert len(result.seqM) == len(result.seqxA) == len(result.seqyA)


@pytest.mark.parametrize(
    "coor_1_path, coor_2_path, chain_1, chain_2, expected_L, expected_rmsd, expected_tm1, expected_tm2",
    [
        (PDB_1JD4, PDB_5M6N, "A", "A", 89, 1.1911595829841628, 0.8425104304046533, 0.875839308266838),
        (PDB_1JD4, PDB_5M6N, "B", "B", 89, 1.408638270685818, 0.9060533678097049, 0.8702330973770631),
        (DOCKQ_MODEL, DOCKQ_NATIVE, "A", "A", 373, 1.120731020605673, 0.9735593880908936, 0.9735593880908936),
        (DOCKQ_MODEL, DOCKQ_NATIVE, "B", "B", 228, 0.9834864903175197, 0.9770609735605671, 0.9770609735605671),
        (PDB_1RXZ_Colabfold, PDB_1RXZ, "B", "A", 245, 0.6864002147683985, 0.9864055216049585, 0.9864055216049585),
        (PDB_1RXZ_Colabfold, PDB_1RXZ, "C", "B", 11, 0.4321917408923342, 0.6843953715215304, 0.6843953715215304),
    ],
)
def test_tmalign_additional_complex_references(
    coor_1_path,
    coor_2_path,
    chain_1,
    chain_2,
    expected_L,
    expected_rmsd,
    expected_tm1,
    expected_tm2,
):
    coor_1 = Coor(coor_1_path)
    coor_2 = Coor(coor_2_path)

    result = tmalign_ca(coor_1, coor_2, chain_1=[chain_1], chain_2=[chain_2], mm=1)

    assert result.L_ali == expected_L
    assert result.rmsd == pytest.approx(expected_rmsd, abs=1e-6)
    assert result.TM1 == pytest.approx(expected_tm1, abs=1e-6)
    assert result.TM2 == pytest.approx(expected_tm2, abs=1e-6)


def test_compute_secondary_structure_single_model_reference():
    coor = Coor(PDB_1Y0M)

    ss_list = TMalign.compute_secondary_structure(coor)
    seqs = coor.get_aa_seq()

    assert len(ss_list) == len(seqs)
    assert len(ss_list[0]["A"]) == len(seqs["A"])
    assert (
        ss_list[0]["A"]
        == "CCEEEEEECCEECCCTTCCEEECTCCEECCCEEECHCCEEEECTTCCCEEECTHCCEEECC"
    )


def test_compute_secondary_structure_multimodel_reference():
    coor = Coor(PDB_1U85)

    ss_list = TMalign.compute_secondary_structure(coor)

    ss_expected = [
        "CCCCECTCECCHTCEECCCTHHHHHHHHHHCCC",
        "CCECCCTCECCHTCEECCCTHHHHHHHHHCCCC",
        "CCCCCCTCECCHTCCECCCTHHHHHHHHHHTCC",
        "CCECCCTCECCHTCCECCCTHHHHHHHHHHTCC",
        "CCCCCCTCECCHTCCECCCTHHHHHHHHHHTCC",
        "CCCCCCTCECCHTCCECCCTHHHHHHHHHCCCC",
        "CCEEECTCEECHTCCECCCTHHHHHHHHHCCCC",
        "CCECCCTCEECHTCEECCCTHHHHHHHHHCCCC",
        "CCCHCCTCECCHTCCECCCTHHHHHHHHHHTCC",
        "CCHCECTCECCHTCCECCCTHHHHHHHHHHTCC",
        "CCCCECTCECCHTCCECCCTHHHHHHHHHCTCC",
        "CCCCECTCEECHTCEECCCTHHHHHHCHHHCCC",
        "CCHHCCTCECCHTCEEECCTHHHHHHHHHCCCC",
        "CCCHCCTCECCHTCCECCCTHHHHHHHHHHHCC",
        "CCCCCCTCECCHTCEEECCTHHHHHHHHHHTCC",
        "CCEEECTCECCHTCEEECCTHHHHHHHHHCCCC",
        "CCCEECTCECCHTCCECCCTHHHHHHHHHHCCC",
        "CCECCCTCECCHTCEECCCTHHHHHHHHHHTCC",
        "CCCCCCTCECCHTCCECCCTHHHHHHHHHHHCC",
        "CCCHCCTCEECHTCEEECCTHHHHHHHHHCHCC",
    ]

    assert len(ss_list) == len(ss_expected)
    for index, expected in enumerate(ss_expected):
        assert ss_list[index]["A"] == expected
