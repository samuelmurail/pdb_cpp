#!/usr/bin/env python3
# coding: utf-8

"""
Tests for _alignement functions
"""

from .datafiles import PDB_1Y0M, PDB_1U85
import pdb_cpp
from pdb_cpp import Coor, TMalign


def test_dssp():
    coor = Coor(PDB_1Y0M)

    SS_list = TMalign.compute_secondary_structure(coor)
    seqs = coor.get_aa_seq()
    print(seqs["A"])
    assert len(SS_list) == len(seqs)

    assert len(SS_list[0]["A"]) == len(seqs["A"])

    assert (
        SS_list[0]["A"]
        == "CCEEEEEECCEECCCTTCCEEECTCCEECCCEEECHCCEEEECTTCCCEEECTHCCEEECC"
    )



def test_multiple_dssp():
    coor = Coor(PDB_1U85)

    SS_list = TMalign.compute_secondary_structure(coor)

    SS_expected = [
        'CCCCECTCECCHTCEECCCTHHHHHHHHHHCCC',
        'CCECCCTCECCHTCEECCCTHHHHHHHHHCCCC',
        'CCCCCCTCECCHTCCECCCTHHHHHHHHHHTCC',
        'CCECCCTCECCHTCCECCCTHHHHHHHHHHTCC',
        'CCCCCCTCECCHTCCECCCTHHHHHHHHHHTCC',
        'CCCCCCTCECCHTCCECCCTHHHHHHHHHCCCC',
        'CCEEECTCEECHTCCECCCTHHHHHHHHHCCCC',
        'CCECCCTCEECHTCEECCCTHHHHHHHHHCCCC',
        'CCCHCCTCECCHTCCECCCTHHHHHHHHHHTCC',
        'CCHCECTCECCHTCCECCCTHHHHHHHHHHTCC',
        'CCCCECTCECCHTCCECCCTHHHHHHHHHCTCC',
        'CCCCECTCEECHTCEECCCTHHHHHHCHHHCCC',
        'CCHHCCTCECCHTCEEECCTHHHHHHHHHCCCC',
        'CCCHCCTCECCHTCCECCCTHHHHHHHHHHHCC',
        'CCCCCCTCECCHTCEEECCTHHHHHHHHHHTCC',
        'CCEEECTCECCHTCEEECCTHHHHHHHHHCCCC',
        'CCCEECTCECCHTCCECCCTHHHHHHHHHHCCC',
        'CCECCCTCECCHTCEECCCTHHHHHHHHHHTCC',
        'CCCCCCTCECCHTCCECCCTHHHHHHHHHHHCC',
        'CCCHCCTCEECHTCEEECCTHHHHHHHHHCHCC',
    ]

    assert len(SS_list) == len(SS_expected)

    for i in range(len(SS_list)):
        assert SS_list[i]["A"] == SS_expected[i]

