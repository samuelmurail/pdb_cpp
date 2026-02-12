#!/usr/bin/env python3
# coding: utf-8


def remove_incomplete_backbone_residues(coor, back_atom=None):
    if back_atom is None:
        back_atom = ["CA", "C", "N", "O"]
    return coor.remove_incomplete_backbone_residues(back_atom)
