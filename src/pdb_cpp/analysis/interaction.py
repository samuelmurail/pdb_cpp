#!/usr/bin/env python3
# coding: utf-8

"""Grouped interaction-analysis helpers.

This module collects contact and interface-oriented utilities under the
analysis namespace.
"""

from __future__ import annotations

from .hbonds import hbonds
from .salt_bridge import salt_bridges
from .sasa import buried_surface_area

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

    This is a semantic alias for :func:`pdb_cpp.analysis.sasa.buried_surface_area`
    exposed under the analysis interaction namespace.
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