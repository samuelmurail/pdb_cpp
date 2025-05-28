#!/usr/bin/env python3
# coding: utf-8

from .core import Alignment, sequence_align
from .data.blosum import BLOSUM62

def align_seq(seq1, seq2, gap_cost=-11, gap_ext=-1, matrix_file='src/pdb_cpp/data/blosum62.txt'):
    """
    Align two sequences using a simple scoring system.

    Parameters
    ----------
    seq_1 : str
        First sequence
    seq_2 : str
        Second sequence
    gap_cost : int, optional
        Cost for opening a gap, by default -11
    gap_ext : int, optional
        Cost for extending a gap, by default 1
    matrix_file : str, optional
        Path to the scoring matrix file, by default 'src/pdb_cpp/data/blosum62.txt'

    Returns
    -------
    tuple
        Tuple containing the aligned sequences (seq1, seq2)
    """

    alignement = sequence_align(
        seq1=seq1,
        seq2=seq2,
        matrix_file='src/pdb_cpp/data/blosum62.txt',
        GAP_COST=gap_cost,
        GAP_EXT=gap_ext)

    return alignement.seq1, alignement.seq2


def print_align_seq(seq_1, seq_2, line_len=80):
    """Print the aligned sequences with a line length of 80 characters.

    Parameters
    ----------
    seq_1 : str
        First sequence
    seq_2 : str
        Second sequence
    line_len : int, optional
        Length of the line, by default 80

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
                # print(seq_1[i], seq_2[i])
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
