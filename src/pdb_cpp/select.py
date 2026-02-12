#!/usr/bin/env python3
# coding: utf-8

__all__ = ["remove_incomplete_backbone_residues"]


def remove_incomplete_backbone_residues(coor, back_atom=None):
    """Remove residues with incomplete backbone atoms.

    Parameters
    ----------
    coor : Coor
        Coordinate object to clean.
    back_atom : list[str], optional
        Backbone atom names to require per residue.

    Returns
    -------
    Coor
        A new Coor object with incomplete residues removed.
    """

    if back_atom is None:
        back_atom = ["CA", "C", "N", "O"]
    return coor.remove_incomplete_backbone_residues(back_atom)
