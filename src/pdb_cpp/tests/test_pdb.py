#!/usr/bin/env python3
# coding: utf-8

"""
Tests for pdb_manip functions
"""

import os
import pytest
import numpy as np
import logging

from pdb_cpp import Coor
from .datafiles import PDB_1Y0M, PDB_2RRI, PDB_5M6N#, PQR_1Y0M, PDB_3FTK


def test_get_pdb(tmp_path):
    """Test get_pdb function."""
    test = Coor(PDB_1Y0M)
    assert test.len == 648
    assert test.model_num == 1

    #assert test.models[0].atom_dict["name_resname_elem"][0, 1] == "THR"
    assert ''.join(test.models[0].resname[0]) == "THR\x00\x00"
    assert test.models[0].resid[0] == 791
    assert test.models[0].uniq_resid[0] == 0
    assert test.models[0].name[0][0] == "N"
    assert test.models[0].num[0] == 1
    assert test.models[0].x[0] == pytest.approx(-1.432, 0.000001)
    assert test.models[0].y[0] == pytest.approx(9.759, 0.000001)
    assert test.models[0].z[0] == pytest.approx(11.436, 0.000001)
    # assert (
    #     test.crystal_pack
    #     == "CRYST1   28.748   30.978   29.753  90.00  92.12  90.00 P 1 21 1      2          \n"
    # )


def test_get_pdb_models(tmp_path):
    """Test get_pdb function."""
    test = Coor(PDB_2RRI)

    assert test.model_num == 20

    assert test.len == 479
    assert test.active_model == 0

    assert test.resname[0] == ['H', 'I', 'S', '\x00', '\x00']
    assert test.models[0].resname[0] == ['H', 'I', 'S', '\x00', '\x00']
    assert test.resid[0] == 1
    assert test.name[0][0] == "N"
    assert test.num[0] == 1
    assert test.x[0] == pytest.approx(-11.432, 0.000001)
    assert test.y[0] == pytest.approx(14.757, 0.000001)
    assert test.z[0] == pytest.approx(-14.63, 0.000001)

    test.active_model = 19
    assert test.active_model == 19

    assert test.resname[0] == ['H', 'I', 'S', '\x00', '\x00']
    assert test.models[19].resid[0] == 1
    assert test.models[19].name[0][0] == "N"
    assert test.name[0][0] == "N"
    assert test.models[19].num[0] == 1
    assert test.models[19].x[0] == pytest.approx(1.574, 0.000001)
    assert test.models[19].y[0] == pytest.approx(17.66, 0.000001)
    assert test.models[19].z[0] == pytest.approx(-0.301, 0.000001)
    # assert (
    #     test.crystal_pack
    #     == "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          \n"
    # )


def test_read_write_pdb(tmp_path, caplog):
    """Test read_file function."""


    test = Coor(PDB_1Y0M)
    assert test.len == 648

    test.write(os.path.join(tmp_path, "test.pdb"))
    captured = caplog.records


    test2 = Coor(os.path.join(tmp_path, "test.pdb"))
    assert test2.len == test.len
    #assert test2.crystal_pack.strip() == test.crystal_pack.strip()

    assert (
        test.models[0].x == test2.models[0].x
        )
    assert (
        test.models[0].resname == test2.models[0].resname
        )


def test_read_write_pdb_models(tmp_path):
    """Test read and write pdb function with several models."""
    test = Coor(PDB_2RRI)

    assert test.model_num == 20

    test.write(os.path.join(tmp_path, "test_2rri.pdb"))

    test_2 = Coor(os.path.join(tmp_path, "test_2rri.pdb"))

    assert test_2.model_num == 20

    assert (
        test.models[0].x == test_2.models[0].x
        )
    assert (
        test.models[0].resname == test_2.models[0].resname
        )


def test_pdb_conect(tmp_path):
    test = Coor(PDB_5M6N)

    assert test.len == 2147
    assert test.model_num == 1

    assert len(test.conect) == 191
    for key in test.conect:
        assert len(test.conect[key]) >= 1

    assert test.conect[369] == [1618]
    assert test.conect[1743] == [1755, 1756, 1782, 1783]
    assert test.conect[1800] == [1759]

    assert test.conect[1629] == [1630, 1666, 1667, 1668]

    test_sel = test.select_atoms("protein")

    assert len(test_sel.conect) == 0

    test_sel_2 = test.select_atoms("resname 7H9")

    assert len(test_sel_2.conect) == 156

    for key in test_sel_2.conect:
        assert len(test_sel_2.conect[key]) >= 1

    assert test_sel_2.conect[1] == [2, 38, 39, 40]
    assert test_sel_2.conect[83] == [82, 84, 122, 123]
    assert test_sel_2.conect[156] == [115]

    test_sel_2.write(os.path.join(tmp_path, "test_5m6n_conect.pdb"))

    test_reload = Coor(os.path.join(tmp_path, "test_5m6n_conect.pdb"))
    assert len(test_reload.conect) == 156
    for key in test_reload.conect:
        assert len(test_reload.conect[key]) >= 1
    assert test_reload.conect[1] == [2, 38, 39, 40]
    assert test_reload.conect[83] == [82, 84, 122, 123]
    assert test_reload.conect[155] == [115]


# def test_get_pdb_bioassembly(tmp_path):
#     """Test get_pdb function."""
#     test = Coor(PDB_3FTK)
#     test.merge_models()

#     assert test.len == 58
#     assert test.model_num == 1

#     assert test.models[0].atom_dict["name_resname_elem"][0, 1] == "ASN"
#     assert test.models[0].resname[0] == "ASN"
#     assert test.models[0].resid[0] == 1
#     assert test.models[0].uniq_resid[0] == 0
#     assert test.models[0].name[0] == "N"
#     assert test.models[0].num[0] == 1
#     assert test.models[0].x[0] == pytest.approx(-8.053, 0.000001)
#     assert test.models[0].y[0] == pytest.approx(2.244, 0.000001)
#     assert test.models[0].z[0] == pytest.approx(10.035, 0.000001)
#     assert (
#         test.models[0].atom_dict["xyz"][0, :]
#         == np.array([-8.053, 2.244, 10.035], dtype=np.float32)
#     ).all()
#     assert (
#         test.crystal_pack
#         == "CRYST1   20.630    4.700   21.009  90.00  92.28  90.00 P 1 21 1      2          \n"
#     )

#     test2 = Coor()
#     test2 = pdb.fetch_BioAssembly("3FTK", index=1)
#     test2.merge_models()
#     test2.compute_chains_CA()

#     assert test2.len == 174
#     assert test2.model_num == 1

#     assert test2.models[0].atom_dict["name_resname_elem"][0, 1] == "ASN"
#     assert test2.models[0].resname[0] == "ASN"
#     assert test2.models[0].resid[0] == 1
#     assert test2.models[0].uniq_resid[0] == 0
#     assert test2.models[0].name[0] == "N"
#     assert test2.models[0].num[0] == 1
#     assert test2.models[0].x[0] == pytest.approx(-8.053, 0.000001)
#     assert test2.models[0].y[0] == pytest.approx(2.244, 0.000001)
#     assert test2.models[0].z[0] == pytest.approx(10.035, 0.000001)

#     assert len(np.unique(test2.models[0].chain)) == 3


# def test_pdb_symmetry_assembly(tmp_path):
#     """Test get_pdb function."""
#     test = Coor(PDB_3FTK)
#     test.merge_models()

#     assert test.len == 58
#     assert test.model_num == 1

#     test = test.add_symmetry()

#     assert test.len == 116
#     assert test.model_num == 1

#     test = test.apply_transformation(index_list=[1, 2])

#     assert test.len == 696
#     assert test.model_num == 1

#     assert len(np.unique(test.chain)) == 1

#     test.compute_chains_CA()

#     assert len(np.unique(test.chain)) == 12

#     test = test.remove_overlap_chain()
#     assert test.len == 431
