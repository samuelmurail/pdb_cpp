#!/usr/bin/env python3
# coding: utf-8

from importlib import resources

from .core import align_chain_permutation as _align_chain_permutation
from .core import sequence_align
from .data.blosum import BLOSUM62

__all__ = ["align_seq", "print_align_seq", "align_chain_permutation"]


def _default_matrix_file():
    return str(resources.files("pdb_cpp.data").joinpath("blosum62.txt"))


def align_seq(seq1, seq2, gap_cost=-11, gap_ext=-1, matrix_file=None):
    """Align two sequences using a simple scoring system.

    Parameters
    ----------
    seq1 : str
        First sequence.
    seq2 : str
        Second sequence.
    gap_cost : int, optional
        Cost for opening a gap.
    gap_ext : int, optional
        Cost for extending a gap.
    matrix_file : str, optional
        Path to the scoring matrix file. If None, uses packaged BLOSUM62.

    Returns
    -------
    tuple
        Tuple containing the aligned sequences (seq1, seq2) and score.
    """
    if matrix_file is None:
        matrix_file = _default_matrix_file()

    alignment = sequence_align(
        seq1=seq1,
        seq2=seq2,
        matrix_file=matrix_file,
        GAP_COST=gap_cost,
        GAP_EXT=gap_ext,
    )

    return alignment.seq1, alignment.seq2, alignment.score


def print_align_seq(seq_1, seq_2, line_len=80):
    """Print the aligned sequences with a fixed line length.

    Parameters
    ----------
    seq_1 : str
        First sequence.
    seq_2 : str
        Second sequence.
    line_len : int, optional
        Length of each output line.

    Returns
    -------
    None
    """

    sim_seq = ""
    for i in range(len(seq_1)):
        if seq_1[i] == seq_2[i]:
            sim_seq += "*"
            continue
        elif seq_1[i] != "-" and seq_2[i] != "-":
            if (seq_1[i], seq_2[i]) in BLOSUM62:
                mut_score = BLOSUM62[seq_1[i], seq_2[i]]
            else:
                mut_score = BLOSUM62[seq_2[i], seq_1[i]]
            if mut_score >= 0:
                sim_seq += "|"
                continue
        sim_seq += " "

    for i in range(1 + len(seq_1) // line_len):
        print(seq_1[i * line_len : (i + 1) * line_len])
        print(sim_seq[i * line_len : (i + 1) * line_len])
        print(seq_2[i * line_len : (i + 1) * line_len])
        print("\n")

    identity = 0
    similarity = 0
    for char in sim_seq:
        if char == "*":
            identity += 1
        if char in ["|", "*"]:
            similarity += 1

    len_1 = len(seq_1.replace("-", ""))
    len_2 = len(seq_2.replace("-", ""))

    print(f"Identity seq1: {identity / len_1 * 100:.2f}%")
    print(f"Identity seq2: {identity / len_2 * 100:.2f}%")

    print(f"Similarity seq1: {similarity / len_1 * 100:.2f}%")
    print(f"Similarity seq2: {similarity / len_2 * 100:.2f}%")

    return


def align_chain_permutation(
    coor_1,
    coor_2,
    back_names=None,
    matrix_file=None,
    frame_ref=0,
):
    """Align structures by permuting chain order and selecting the best RMSD.

    Parameters
    ----------
    coor_1 : Coor
        First coordinate object.
    coor_2 : Coor
        Second coordinate object.
    back_names : list[str], optional
        Backbone atom names to use.
    matrix_file : str, optional
        Path to the scoring matrix file. If None, uses packaged BLOSUM62.
    frame_ref : int, optional
        Reference frame index in coor_2.

    Returns
    -------
    tuple
        RMSD list and index mappings from the best permutation.
    """
    if back_names is None:
        back_names = ["C", "N", "O", "CA"]

    if matrix_file is None:
        matrix_file = _default_matrix_file()

    return _align_chain_permutation(
        coor_1,
        coor_2,
        back_names,
        matrix_file,
        frame_ref,
    )
