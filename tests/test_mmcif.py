#!/usr/bin/env python3
# coding: utf-8

"""Tests for mmCIF parsing and writing."""

import os
import pytest

from pdb_cpp import Coor
from .datafiles import MMCIF_1Y0M, MMCIF_2RRI, MMCIF_9X0F


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
