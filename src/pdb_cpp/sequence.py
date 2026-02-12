#!/usr/bin/env python3
# coding: utf-8

from .core import Coor


def get_aa_seq(self, gap_in_seq=True, frame=0):
    """Return the amino acid sequence of the selection.

    Parameters
    ----------
    gap_in_seq : bool, optional
        Whether to insert gaps for missing residues.
    frame : int, optional
        Frame index to use for the selection.

    Returns
    -------
    dict
        Mapping of chain ID to sequence.
    """
    sequences = self.get_aa_sequences(gap_in_seq=gap_in_seq, frame=frame)
    if not sequences:
        return {}
    ca_index = self.get_index_select("name CA")
    chains = self.chain
    seq_dict = {}

    chain_keys = []
    seen = set()
    for idx in ca_index:
        chain_id = ""
        for letter in chains[idx]:
            if letter != "\x00" and letter != " ":
                chain_id += letter
        if chain_id not in seen:
            chain_keys.append(chain_id)
            seen.add(chain_id)

    for chain_id, seq in zip(chain_keys, sequences):
        seq_dict[chain_id] = seq
    return seq_dict


def get_aa_DL_seq(self, gap_in_seq=True, frame=0):
    """Return the amino acid sequence with D-residues in lowercase.

    Parameters
    ----------
    gap_in_seq : bool, optional
        Whether to insert gaps for missing residues.
    frame : int, optional
        Frame index to use for the selection.

    Returns
    -------
    dict
        Mapping of chain ID to sequence.
    """
    sequences = self.get_aa_sequences_dl(gap_in_seq=gap_in_seq, frame=frame)
    if not sequences:
        return {}
    ca_index = self.get_index_select("name CA")
    chains = self.chain
    seq_dict = {}

    chain_keys = []
    seen = set()
    for idx in ca_index:
        chain_id = ""
        for letter in chains[idx]:
            if letter != "\x00" and letter != " ":
                chain_id += letter
        if chain_id not in seen:
            chain_keys.append(chain_id)
            seen.add(chain_id)

    for chain_id, seq in zip(chain_keys, sequences):
        seq_dict[chain_id] = seq
    return seq_dict


Coor.get_aa_seq = get_aa_seq
Coor.get_aa_DL_seq = get_aa_DL_seq
