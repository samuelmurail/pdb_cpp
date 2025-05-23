#!/usr/bin/env python3
# coding: utf-8

"""
Tests for pdb_manip functions
"""

import os
import pytest
import numpy as np

from pdb_cpp import Coor
# from pdb_cpp import select as select

from .datafiles import PDB_1Y0M, PDB_2RRI


# def test_parse_selection():
#     """Test parse_selection function."""

#     selection = "name CA and resname ALA"
#     selec = select.parse_selection(selection)
#     assert selec == [["name", "CA"], "and", ["resname", "ALA"]]

#     selection = "(name CA and resname ALA) or (name C and resname GLY)"
#     selec = select.parse_selection(selection)
#     assert selec == [
#         [["name", "CA"], "and", ["resname", "ALA"]],
#         "or",
#         [["name", "C"], "and", ["resname", "GLY"]],
#     ]

#     selection = "name CA and within 8.0 of not resname GLY TRP"
#     selec = select.parse_selection(selection)
#     assert selec == [
#         ["name", "CA"],
#         "and",
#         ["within", "8.0", "of", ["not", ["resname", "GLY", "TRP"]]],
#     ]

#     selection = "name CA and x < 10 and y >= -10"
#     selec = select.parse_selection(selection)
#     assert selec == [["name", "CA"], "and", ["x", "<", "10"], "and", ["y", ">=", "-10"]]

#     selection = "chain A and protein"
#     selec = select.parse_selection(selection)
#     assert selec == [
#         ["chain", "A"],
#         "and",
#         [
#             "resname",
#             "GLY",
#             "HIS",
#             "HSP",
#             "HSE",
#             "HSD",
#             "HIP",
#             "HIE",
#             "HID",
#             "ARG",
#             "LYS",
#             "ASP",
#             "ASPP",
#             "ASH",
#             "GLU",
#             "GLUP",
#             "GLH",
#             "SER",
#             "THR",
#             "ASN",
#             "GLN",
#             "CYS",
#             "SEC",
#             "PRO",
#             "ALA",
#             "ILE",
#             "PHE",
#             "TYR",
#             "TRP",
#             "VAL",
#             "LEU",
#             "MET",
#             "DAL",
#             "DAR",
#             "DSG",
#             "DAS",
#             "DCY",
#             "DGN",
#             "DGL",
#             "DHI",
#             "DIL",
#             "DLE",
#             "DLY",
#             "DME",
#             "MED",
#             "DPH",
#             "DPN",
#             "DPR",
#             "DSE",
#             "DSN",
#             "DTH",
#             "DTR",
#             "DTY",
#             "DVA",
#         ],
#     ]

#     selection = "chain A and backbone"
#     selec = select.parse_selection(selection)
#     assert selec == [
#         ["chain", "A"],
#         "and",
#         [
#             "resname",
#             "GLY",
#             "HIS",
#             "HSP",
#             "HSE",
#             "HSD",
#             "HIP",
#             "HIE",
#             "HID",
#             "ARG",
#             "LYS",
#             "ASP",
#             "ASPP",
#             "ASH",
#             "GLU",
#             "GLUP",
#             "GLH",
#             "SER",
#             "THR",
#             "ASN",
#             "GLN",
#             "CYS",
#             "SEC",
#             "PRO",
#             "ALA",
#             "ILE",
#             "PHE",
#             "TYR",
#             "TRP",
#             "VAL",
#             "LEU",
#             "MET",
#             "DAL",
#             "DAR",
#             "DSG",
#             "DAS",
#             "DCY",
#             "DGN",
#             "DGL",
#             "DHI",
#             "DIL",
#             "DLE",
#             "DLY",
#             "DME",
#             "MED",
#             "DPH",
#             "DPN",
#             "DPR",
#             "DSE",
#             "DSN",
#             "DTH",
#             "DTR",
#             "DTY",
#             "DVA",
#         ],
#         "and",
#         ["name", "N", "CA", "C", "O"],
#     ]

#     selection = "chain A and noh"
#     selec = select.parse_selection(selection)
#     assert selec == [["chain", "A"], "and", "not", ["name", "H*"]]


def test_select_atoms_model():
    """Test select_atoms function."""
    test = Coor(PDB_1Y0M)
    model = test.get_Models(0)

    assert model.size() == 648
    selec = "name CA and resname ALA"
    new = model.select_atoms(selec)
    assert  sum(new) == 4 

    selec = "backbone and resid >= 796 and resid <= 848"
    new = model.select_atoms(selec)
    assert  sum(new) == 214

    selec = "protein and resid >= 796 and resid <= 848"
    new = model.select_atoms(selec)
    assert  sum(new) == 463

    selec = "resname HOH and chain A"
    new = model.select_atoms(selec)
    assert  sum(new) == 122

    selec = "not resname HOH and chain A"
    new = model.select_atoms(selec)
    print(sum(new))
    assert  sum(new) == 526

    selec = "name CA and chain A and resid >= 796 and resid <= 848"
    new = model.select_atoms(selec)
    assert  sum(new) == 53

    selec = "name CA and chain A and residue >= 6 and residue <= 58"
    new = model.select_atoms(selec)
    assert  sum(new) == 53

    selec = "name CA CX CY and chain A B C and residue >= 6 and residue <= 58"
    new = model.select_atoms(selec)
    assert  sum(new) == 53

    selec = "name CA and chain A and residue >= 6 and residue <= 58 and residue != 30"
    new = model.select_atoms(selec)
    assert  sum(new) == 52

    selec = "resname HOH and chain A and x > 20"
    new = model.select_atoms(selec)
    assert  sum(new) == 56

    selec = "resname HOH and chain A and x >= 20.0"
    new = model.select_atoms(selec)
    assert  sum(new) ==  56

    selec = "name N and residue == 0"
    new = model.select_atoms(selec)
    assert  sum(new) ==  1

    selec = "name N CA and residue 0 1 2"
    new = model.select_atoms(selec)
    assert  sum(new) == 6


def test_select_atoms_multi_frame():
    """Test select_atoms function."""
    test = Coor(PDB_2RRI)
    assert test.len == 479

    

    selec = "name N CA and residue > 20 and residue < 80"
    new = test.select_atoms(selec)
    new_model_10 = new.get_Models(10)
    assert new.len == 16
    assert new_model_10.len == 16

    selec = "name N CA and residue > -20 and residue < 80"
    new = test.select_atoms(selec)
    new_model_10 = new.get_Models(10)
    assert new.len == 58
    assert new_model_10.len == 58

    selec = "name N CA and residue > 20 and residue < 80"
    new = test.select_atoms(selec, frame=19)
    new_model_10 = new.get_Models(10)
    assert new.len == 16
    assert new_model_10.len == 16

    selec = "x > 10"
    new = test.select_atoms(selec)
    new_model_10 = new.get_Models(10)
    assert new.len == 57
    assert new_model_10.len == 57

    selec = "x > 10"
    new = test.select_atoms(selec, frame=10)
    new_model_10 = new.get_Models(10)
    assert new.len == 58
    assert new_model_10.len == 58

    selec = "noh"
    new = test.select_atoms(selec)
    new_model_10 = new.get_Models(10)
    assert new.len == 237
    assert new_model_10.len == 237


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
