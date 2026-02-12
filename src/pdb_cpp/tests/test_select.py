#!/usr/bin/env python3
# coding: utf-8

"""
Tests for pdb_cpp functions
"""

import os
import pytest
import numpy as np

from pdb_cpp import Coor

# from pdb_cpp import select as select

from .datafiles import PDB_1U85, PDB_1Y0M, PDB_2RRI


def test_select_atoms():
    """Test select_atoms function."""
    test = Coor(PDB_1Y0M)
    assert test.len == 648

    selec = "name CA and resname ALA"
    new = test.select_atoms(selec)
    assert new.len == 4

    selec = "backbone and resid >= 796 and resid <= 848"
    new = test.select_atoms(selec)
    assert new.len == 214

    selec = "protein and resid >= 796 and resid <= 848"
    new = test.select_atoms(selec)
    assert new.len == 463

    selec = "resname HOH and chain A"
    new = test.select_atoms(selec)
    assert new.len == 122

    selec = "name CA and chain A and resid >= 796 and resid <= 848"
    new = test.select_atoms(selec)
    assert new.len == 53

    selec = "name CA and chain A and residue >= 6 and residue <= 58"
    new = test.select_atoms(selec)
    assert new.len == 53

    selec = "name CA CX CY and chain A B C and residue >= 6 and residue <= 58"
    new = test.select_atoms(selec)
    assert new.len == 53

    selec = "name CA and chain A and residue >= 6 and residue <= 58 and residue != 30"
    new = test.select_atoms(selec)
    assert new.len == 52

    selec = "resname HOH and chain A and x > 20"
    new = test.select_atoms(selec)
    assert new.len == 56

    selec = "resname HOH and chain A and x >= 20.0"
    new = test.select_atoms(selec)
    assert new.len == 56

    selec = "name N and residue == 0"
    print(test.models[0].residue)
    print(test.models[0].name)
    new = test.select_atoms(selec)
    assert new.len == 1

    selec = "name N CA and residue 0 1 2"
    new = test.select_atoms(selec)
    assert new.len == 6


def test_select_atoms_multi_frame():
    """Test select_atoms function."""
    test = Coor(PDB_2RRI)
    assert test.len == 479

    selec = "name N CA and residue > 20 and residue < 80"
    new = test.select_atoms(selec)
    assert new.len == 16
    assert new.models[10].len == 16

    selec = "name N CA and residue > -20 and residue < 80"
    new = test.select_atoms(selec)
    assert new.len == 58
    assert new.models[10].len == 58

    selec = "name N CA and residue > 20 and residue < 80"
    new = test.select_atoms(selec, frame=19)
    assert new.len == 16
    assert new.models[10].len == 16

    selec = "x > 10"
    new = test.select_atoms(selec)
    assert new.len == 57
    assert new.models[10].len == 57

    selec = "x > 10"
    new = test.select_atoms(selec, frame=10)
    assert new.len == 58
    assert new.models[10].len == 58

    selec = "noh"
    new = test.select_atoms(selec)
    assert new.len == 237
    assert new.models[10].len == 237


def test_select_atoms_within(tmp_path):
    """Test select_atoms function with within selection."""

    test = Coor(PDB_1Y0M)

    selec = "name CA and chain A"
    new = test.select_atoms(selec)
    assert new.len == 61

    selec = "name CA and within 5.0 of resname HOH and chain A"
    new = test.select_atoms(selec)
    assert new.len == 48

    selec = "name CA and within 1.0 of resname HOH and chain A"
    new = test.select_atoms(selec)
    assert new.len == 0

    selec = "name CX and within 10.0 of resname HOH and chain A"
    new = test.select_atoms(selec)
    assert new.len == 0

    selec = "name CA and not within 5 of resname HOH and chain A"
    new = test.select_atoms(selec)
    assert new.len == 13

    selec = "name CA and within 5 of not resname HOH ALA GLY SER TRP THR PRO TYR PHE GLU ASP HIS ARG LYS"
    new = test.select_atoms(selec)
    assert new.len == 49


def test_select_atoms_within_multi_frame():
    """Test select_atoms function."""
    test = Coor(PDB_2RRI)
    assert test.len == 479

    selec = "residue 20"
    new = test.select_atoms(selec)
    assert new.len == 22

    selec = "within 5 of residue 20"
    new = test.select_atoms(selec)
    assert new.len == 96
    assert new.models[19].len == 96

    selec = "within 5 of residue 20"
    new = test.select_atoms(selec, frame=15)
    assert new.len == 99
    assert new.models[10].len == 99


def test_select_atoms_conect_update():
    """Ensure CONECT records are renumbered after selection."""
    test = Coor(PDB_1U85)

    assert 529 in test.conect
    assert 139 in test.conect

    new = test.select_atoms("num 139 529")
    assert new.len == 2

    assert new.num == [1, 2]
    assert new.conect[1] == [2]
    assert new.conect[2] == [1]
