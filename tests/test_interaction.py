#!/usr/bin/env python3
# coding: utf-8

from pdb_cpp import Coor
from pdb_cpp.analysis import hbonds as analysis_hbonds
from pdb_cpp.analysis import interaction as analysis_interaction
from pdb_cpp.analysis import salt_bridge as analysis_salt_bridge

from .datafiles import CIF_1A0A, MMCIF_1Y0M, PDB_1A2K_MODEL


def test_interaction_hbonds_matches_analysis_hbonds_module():
    coor = Coor(MMCIF_1Y0M)
    ref = analysis_hbonds.hbonds(coor, angle_cutoff=120.0)
    grouped = analysis_interaction.hbonds(coor, angle_cutoff=120.0)
    assert len(grouped[0]) == len(ref[0])


def test_interaction_salt_bridges_matches_analysis_salt_bridge_module():
    coor = Coor(CIF_1A0A)
    ref = analysis_salt_bridge.salt_bridges(coor, cation_sel="protein", anion_sel="nucleic")
    grouped = analysis_interaction.salt_bridges(coor, cation_sel="protein", anion_sel="nucleic")
    assert len(grouped[0]) == len(ref[0])


def test_interaction_interface_sasa_matches_analysis_wrapper():
    coor = Coor(PDB_1A2K_MODEL)
    result = analysis_interaction.interface_sasa(
        coor,
        receptor_sel="chain A",
        ligand_sel="chain B",
    )
    assert len(result) == coor.model_num
    assert "buried_surface" in result[0]
    assert "interface_area" in result[0]