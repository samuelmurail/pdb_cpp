#!/usr/bin/env python3
# coding: utf-8

"""Grouped interaction-analysis helpers.

This module collects contact and interface-oriented utilities that were
previously spread across separate modules.
"""

from __future__ import annotations

from .analysis import buried_surface_area
from .hbond import hbonds
from .salt_bridge import salt_bridges

__all__ = ["hbonds", "salt_bridges", "interface_sasa"]


def interface_sasa(
    coor,
    receptor_sel="protein",
    ligand_sel="not protein",
    probe_radius=1.4,
    n_points=960,
    include_hydrogen=False,
    by_residue=False,
):
    """Compute interface SASA / buried surface area for two selections.

    This is a semantic alias for :func:`pdb_cpp.analysis.buried_surface_area`
    exposed under the interaction namespace.
    """
    return buried_surface_area(
        coor,
        receptor_sel=receptor_sel,
        ligand_sel=ligand_sel,
        probe_radius=probe_radius,
        n_points=n_points,
        include_hydrogen=include_hydrogen,
        by_residue=by_residue,
    )