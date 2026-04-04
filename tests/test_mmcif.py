#!/usr/bin/env python3
# coding: utf-8

"""Tests for mmCIF parsing and writing."""

import os
import pytest

from pdb_cpp import Coor
from .datafiles import MMCIF_1Y0M, MMCIF_2RRI, MMCIF_9X0F, PDB_5M6N, PDB_1U85, CIF_1A0A, CIF_FOLD_2026_DNAPROT_MODEL


def test_read_mmcif_basic():
    test = Coor(MMCIF_1Y0M)
    assert test.len == 648
    assert test.model_num == 1

    assert "".join(test.models[0].resname[0]) == "THR\x00\x00"
    assert test.models[0].resid[0] == 791
    assert test.models[0].uniq_resid[0] == 0
    assert test.models[0].name[0][0] == "N"
    assert test.models[0].num[0] == 1
    assert test.models[0].x[0] == pytest.approx(-1.432, 0.000001)
    assert test.models[0].y[0] == pytest.approx(9.759, 0.000001)
    assert test.models[0].z[0] == pytest.approx(11.436, 0.000001)


def test_read_mmcif_write_pdb(tmp_path):
    test = Coor(MMCIF_1Y0M)
    out_path = os.path.join(tmp_path, "test_mmcif.pdb")
    test.write(out_path)

    test_2 = Coor(out_path)
    assert test_2.len == test.len
    assert test_2.model_num == test.model_num
    assert test_2.models[0].resname == test.models[0].resname
    assert test_2.models[0].resid[0] == test.models[0].resid[0]
    assert test_2.models[0].num[0] == test.models[0].num[0]
    assert test_2.models[0].x[0] == pytest.approx(test.models[0].x[0], 0.000001)
    assert test_2.models[0].y[0] == pytest.approx(test.models[0].y[0], 0.000001)
    assert test_2.models[0].z[0] == pytest.approx(test.models[0].z[0], 0.000001)


def test_read_mmcif_write_mmcif(tmp_path):
    test = Coor(MMCIF_1Y0M)
    out_path = os.path.join(tmp_path, "test_mmcif.cif")
    test.write(out_path)

    test_2 = Coor(out_path)
    assert test_2.len == test.len
    assert test_2.model_num == test.model_num
    assert test_2.models[0].resid[0] == test.models[0].resid[0]
    assert test_2.models[0].num[0] == test.models[0].num[0]
    assert test_2.models[0].resname == test.models[0].resname


def test_read_mmcif_write_pdb_models(tmp_path):
    test = Coor(MMCIF_2RRI)

    assert test.model_num == 20
    assert test.len == 479

    out_path = os.path.join(tmp_path, "test_2rri.pdb")
    test.write(out_path)

    test_2 = Coor(out_path)
    assert test_2.model_num == 20
    assert test_2.len == 479

    assert test_2.models[0].resid[0] == test.models[0].resid[0]
    assert test_2.models[0].num[0] == test.models[0].num[0]


def test_read_mmcif_write_mmcif_models(tmp_path):
    test = Coor(MMCIF_2RRI)

    assert test.model_num == 20
    assert test.len == 479

    out_path = os.path.join(tmp_path, "test_2rri.cif")
    test.write(out_path)

    test_2 = Coor(out_path)
    assert test_2.model_num == 20
    assert test_2.len == 479

    assert test_2.models[0].resid[0] == test.models[0].resid[0]


def test_read_mmcif_9x0f():
    test = Coor(MMCIF_9X0F)

    assert test.model_num == 1
    assert test.len == 42552

    selected = test.select_atoms("within 10.0 of chain A")
    assert selected.len > 0


def test_mmcif_conect_roundtrip(tmp_path):
    """PDB CONECT records survive a PDB -> mmCIF -> PDB roundtrip."""
    coor = Coor(PDB_5M6N)
    assert len(coor.conect) == 191
    assert coor.conect[369] == [1618]
    assert sorted(coor.conect[1743]) == [1755, 1756, 1782, 1783]

    # Write to mmCIF
    cif_path = os.path.join(tmp_path, "5m6n.cif")
    coor.write(cif_path)

    # Read back mmCIF
    coor2 = Coor(cif_path)
    assert len(coor2.conect) == 191

    # Check all bonds match
    for k, v in coor.conect.items():
        assert sorted(coor2.conect[k]) == sorted(v), f"Mismatch at atom {k}"


def test_mmcif_conect_small(tmp_path):
    """CONECT roundtrip for a small structure (1u85 with zinc coordination)."""
    coor = Coor(PDB_1U85)
    assert 529 in coor.conect

    cif_path = os.path.join(tmp_path, "1u85.cif")
    coor.write(cif_path)

    coor2 = Coor(cif_path)
    assert len(coor2.conect) == len(coor.conect)

    for k, v in coor.conect.items():
        assert sorted(coor2.conect[k]) == sorted(v)


def test_mmcif_conect_selection(tmp_path):
    """CONECT records are renumbered correctly after selection + mmCIF roundtrip."""
    coor = Coor(PDB_1U85)
    sel = coor.select_atoms("num 139 529")
    assert sel.conect[1] == [2]
    assert sel.conect[2] == [1]

    # Write selection to mmCIF, read back
    cif_path = os.path.join(tmp_path, "sel.cif")
    sel.write(cif_path)

    sel2 = Coor(cif_path)
    assert sel2.conect[1] == [2]
    assert sel2.conect[2] == [1]


def test_mmcif_no_conect():
    """mmCIF without _struct_conn should have empty conect."""
    coor = Coor(MMCIF_1Y0M)
    assert len(coor.conect) == 0


def test_read_mmcif_protein_dna_complex_1a0a():
    """1A0A.cif is a protein-DNA complex: two protein and two DNA chains plus water.

    label_asym_id: A,B = DNA; C,D = protein; E-H = crystallographic water.
    """
    coor = Coor(CIF_1A0A)
    assert coor.len == 1767
    assert coor.model_num == 1

    # Protein chains C and D are returned by get_aa_seq()
    prot_chains = set(coor.get_aa_seq().keys())
    assert prot_chains == {"C", "D"}

    # All four macromolecular chains are returned by get_aa_na_seq()
    macro_chains = set(coor.get_aa_na_seq().keys())
    assert macro_chains == {"A", "B", "C", "D"}

    # Atom-level selections
    assert coor.select_atoms("protein").len == 996
    assert coor.select_atoms("nucleic").len == 691


def test_read_mmcif_protein_dna_model_fold2026():
    """fold_2026_03_10_11_53_model_4.cif is an AlphaFold3 model of a protein-DNA complex.

    label_asym_id: A,B = protein; C,D = DNA.
    """
    coor = Coor(CIF_FOLD_2026_DNAPROT_MODEL)
    assert coor.len == 1695
    assert coor.model_num == 1

    prot_chains = set(coor.get_aa_seq().keys())
    assert prot_chains == {"A", "B"}

    macro_chains = set(coor.get_aa_na_seq().keys())
    assert macro_chains == {"A", "B", "C", "D"}

    assert coor.select_atoms("protein").len == 996
    assert coor.select_atoms("nucleic").len == 699
