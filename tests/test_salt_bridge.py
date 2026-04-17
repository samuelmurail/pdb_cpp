#!/usr/bin/env python3
# coding: utf-8

from pdb_cpp import Coor
from pdb_cpp.analysis import salt_bridge as analysis_salt_bridge

from .datafiles import CIF_1A0A, MMCIF_1Y0M


def _sb_key(sb):
    return (
        sb.cation_chain,
        sb.cation_resid,
        sb.cation_resname,
        sb.cation_name,
        sb.anion_chain,
        sb.anion_resid,
        sb.anion_resname,
        sb.anion_name,
    )


def test_salt_bridges_returns_list_per_model():
    coor = Coor(MMCIF_1Y0M)
    result = analysis_salt_bridge.salt_bridges(coor)
    assert isinstance(result, list)
    assert len(result) == coor.model_num


def test_analysis_salt_bridge_module_imports_work():
    coor = Coor(MMCIF_1Y0M)
    result = analysis_salt_bridge.salt_bridges(coor)
    assert isinstance(result, list)
    assert len(result) == coor.model_num


def test_salt_bridges_protein_protein_1y0m():
    coor = Coor(MMCIF_1Y0M)
    sb_list = analysis_salt_bridge.salt_bridges(coor, cation_sel="protein", anion_sel="protein")[0]
    assert len(sb_list) > 0
    assert all(sb.distance <= 4.0 for sb in sb_list)


def test_salt_bridges_protein_to_nucleic_1a0a():
    coor = Coor(CIF_1A0A)
    sb_list = analysis_salt_bridge.salt_bridges(coor, cation_sel="protein", anion_sel="nucleic")[0]
    keys = {_sb_key(sb) for sb in sb_list}

    assert len(keys) == 20
    assert ("C", 46, "ARG", "NH2", "A", 4, "DC", "OP2") in keys
    assert ("D", 112, "ARG", "NH2", "B", 23, "DC", "OP1") in keys


def test_salt_bridges_nucleic_to_nucleic_canonical_empty_1a0a():
    coor = Coor(CIF_1A0A)
    sb_list = analysis_salt_bridge.salt_bridges(coor, cation_sel="nucleic", anion_sel="nucleic")[0]
    assert sb_list == []