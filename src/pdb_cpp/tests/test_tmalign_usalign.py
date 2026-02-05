#!/usr/bin/env python3
# coding: utf-8

"""Tests for US-align based TM-align integration."""

import pytest

from .datafiles import PDB_1Y0M, PDB_2RRI
from pdb_cpp import Coor
from pdb_cpp.core import tmalign_ca


def test_tmalign_self_alignment():
    coor = Coor(PDB_1Y0M)

    result = tmalign_ca(coor, coor, chain_1=["A"], chain_2=["A"])

    assert 0.95 <= result.TM1 <= 1.0
    assert 0.95 <= result.TM2 <= 1.0
    assert 0.0 <= result.TM_ali <= 1.0
    assert result.rmsd >= 0.0
    assert result.L_ali >= 0
    assert len(result.seqM) == len(result.seqxA) == len(result.seqyA)


def test_tmalign_basic_output_ranges():
    coor_1 = Coor(PDB_1Y0M)
    coor_2 = Coor(PDB_2RRI)

    result = tmalign_ca(coor_1, coor_2, chain_1=["A"], chain_2=["A"])

    assert 0.0 <= result.TM1 <= 1.0
    assert 0.0 <= result.TM2 <= 1.0
    assert 0.0 <= result.TM_ali <= 1.0
    assert result.rmsd >= 0.0
    assert result.L_ali >= 0
    assert len(result.seqM) == len(result.seqxA) == len(result.seqyA)
