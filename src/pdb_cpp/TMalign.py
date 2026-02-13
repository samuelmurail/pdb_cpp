#!/usr/bin/env python3
# coding: utf-8


from .core import compute_SS

__all__ = ["compute_secondary_structure"]


def compute_secondary_structure(coor, **kwargs):
    """Compute secondary structure using the TM-align core.

    You should cite the original TM-align paper if you use this function:
    Zhang, Y., & Skolnick, J. NAR (2005).

    Parameters
    ----------
    coor : Coor
        The Coor object containing the protein structure.
    **kwargs : dict
        Additional keyword arguments for the compute_SS function.

    Returns
    -------
    list
        List of per-chain secondary-structure dictionaries per model.
    """
    ss_list = compute_SS(coor, **kwargs)
    uniq_chains = coor.get_uniq_chain()

    SS_new_list = []

    for ss_frame in ss_list:
        ss_dict = {}
        # print(f"Secondary structure: {ss_sec}")

        for chain, seq in zip(uniq_chains, ss_frame):
            # print(f"Chain: {chain}, Sequence: {seq}")
            new_chain = ""
            for letter in chain:
                if letter != "\x00":
                    new_chain += letter
            ss_dict[new_chain] = seq
        SS_new_list.append(ss_dict)

    return SS_new_list
