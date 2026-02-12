#!/usr/bin/env python3
# coding: utf-8

"""Thin wrappers around core sequence helpers."""

from .core import Coor


def get_aa_seq(coor, gap_in_seq=True, frame=0):
    """Return the amino acid sequence of the selection.

    Parameters
    ----------
    coor : Coor
        Coordinate object.
    gap_in_seq : bool, optional
        Whether to insert gaps for missing residues.
    frame : int, optional
        Frame index to use for the selection.

    Returns
    -------
    dict
        Mapping of chain ID to sequence.
    """
    return coor.get_aa_seq(gap_in_seq=gap_in_seq, frame=frame)


def get_aa_DL_seq(coor, gap_in_seq=True, frame=0):
    """Return the amino acid sequence with D-residues in lowercase.

    Parameters
    ----------
    coor : Coor
        Coordinate object.
    gap_in_seq : bool, optional
        Whether to insert gaps for missing residues.
    frame : int, optional
        Frame index to use for the selection.

    Returns
    -------
    dict
        Mapping of chain ID to sequence.
    """
    return coor.get_aa_DL_seq(gap_in_seq=gap_in_seq, frame=frame)


__all__ = ["get_aa_seq", "get_aa_DL_seq", "Coor"]
