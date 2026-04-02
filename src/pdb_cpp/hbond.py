#!/usr/bin/env python3
# coding: utf-8

"""Hydrogen-bond computation using Baker & Hubbard geometric criteria.

Virtual hydrogen positions are reconstructed from heavy-atom geometry when H
atoms are absent from the coordinate file.  Backbone N-H positions are computed
from the preceding C and the current CA/N atoms; sidechain hydrogens are placed
along the extension of the bond from the nearest heavy-atom neighbour.

Criteria (all must be satisfied):
    D···A distance  < dist_DA_cutoff  (default 3.5 Å)
    H···A distance  < dist_HA_cutoff  (default 2.5 Å)
    D−H···A angle   > angle_cutoff    (default 90°)

Reference
---------
Baker EN & Hubbard RE (1984) Hydrogen bonding in globular proteins.
*Prog Biophys Mol Biol* **44** 97-179.
"""

import logging

from .core import compute_hbonds as _compute_hbonds

__all__ = ["hbonds"]

logger = logging.getLogger(__name__)


def hbonds(
    coor,
    donor_sel="protein",
    acceptor_sel="protein",
    dist_DA_cutoff=3.5,
    dist_HA_cutoff=2.5,
    angle_cutoff=90.0,
):
    """Compute hydrogen bonds between two selections for every model in *coor*.

    Parameters
    ----------
    coor : Coor
        Coordinate object (one or more models / frames).
    donor_sel : str, optional
        Atom selection string for donor atoms (default: ``"protein"``).
        Use ``"protein or nucleic"`` or ``"nucleic"`` to include nucleic acids.
    acceptor_sel : str, optional
        Atom selection string for acceptor atoms (default: ``"protein"``).
    dist_DA_cutoff : float, optional
        Maximum donor-heavy to acceptor distance in Å (default 3.5).
    dist_HA_cutoff : float, optional
        Maximum hydrogen to acceptor distance in Å (default 2.5).
    angle_cutoff : float, optional
        Minimum D−H···A angle in degrees (default 90).

    Returns
    -------
    list[list[HBond]]
        One list of :class:`~pdb_cpp.core.HBond` objects per model frame.
        Each :class:`~pdb_cpp.core.HBond` has the following read-only attributes:

        * ``donor_resid``       – unique residue ID of the donor
        * ``donor_resname``     – residue name of the donor
        * ``donor_chain``       – chain ID of the donor
        * ``donor_heavy_name``  – heavy donor atom name (e.g. ``"N"``, ``"OG"``)
        * ``donor_h_name``      – hydrogen atom name (actual or virtual)
        * ``donor_heavy_xyz``   – (x, y, z) of the donor heavy atom
        * ``donor_h_xyz``       – (x, y, z) of the H (actual or reconstructed)
        * ``acceptor_resid``    – unique residue ID of the acceptor
        * ``acceptor_resname``  – residue name of the acceptor
        * ``acceptor_chain``    – chain ID of the acceptor
        * ``acceptor_name``     – acceptor atom name (e.g. ``"O"``, ``"OD1"``)
        * ``acceptor_xyz``      – (x, y, z) of the acceptor atom
        * ``dist_DA``           – D···A distance (Å)
        * ``dist_HA``           – H···A distance (Å)
        * ``angle_DHA``         – D−H···A angle (degrees)
    """
    results = []
    for frame_idx, model in enumerate(coor.models):
        donor_model    = coor.select_atoms(donor_sel,    frame_idx).models[frame_idx]
        acceptor_model = coor.select_atoms(acceptor_sel, frame_idx).models[frame_idx]
        hb_list = _compute_hbonds(
            donor_model,
            acceptor_model,
            model,
            dist_DA_cutoff,
            dist_HA_cutoff,
            angle_cutoff,
        )
        results.append(hb_list)
    return results
