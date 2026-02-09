#!/usr/bin/env python3
# coding: utf-8

from .core import Coor


def get_aa_seq(self, gap_in_seq=True, frame=0):
    """Return the amino acid sequence of the selection."""
    uniq_chains = self.get_uniq_chain()
    sequences = self.get_aa_sequences(gap_in_seq=gap_in_seq, frame=frame)
    seq_dict = {}

    for chain, seq in zip(uniq_chains, sequences):
        new_chain = ""
        for letter in chain:
            if letter != "\x00":
                new_chain += letter
        seq_dict[new_chain] = seq
    return seq_dict


def get_aa_DL_seq(self, gap_in_seq=True, frame=0):
    """Return the amino acid sequence with D-residues in lowercase."""
    uniq_chains = self.get_uniq_chain()
    sequences = self.get_aa_sequences_dl(gap_in_seq=gap_in_seq, frame=frame)
    seq_dict = {}

    for chain, seq in zip(uniq_chains, sequences):
        new_chain = ""
        for letter in chain:
            if letter != "\x00":
                new_chain += letter
        seq_dict[new_chain] = seq
    return seq_dict


Coor.get_aa_seq = get_aa_seq
Coor.get_aa_DL_seq = get_aa_DL_seq
