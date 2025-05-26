#!/usr/bin/env python3
# coding: utf-8



from .core import compute_SS

def compute_secondary_structure(coor, **kwargs):
    """
    Compute the secondary structure of a protein structure.

    This function uses the TMalign algorithm to compute the secondary structure elements

    You should cite the original TMalign paper if you use this function in your work:
    Zhang, Y., & Skolnick, J. NAR (2005). TM-align: a protein structure alignment
    algorithm based on the TM-score.
    https://doi.org/10.1093/nar/gki524

    Parameters
    ----------
    coor : Coor
        The Coor object containing the protein structure.
    **kwargs : dict
        Additional keyword arguments for the compute_SS function.

    Returns
    -------
    list
        A list of secondary structure elements.
    """
    ss_list = compute_SS(coor, **kwargs)
    uniq_chains = coor.get_uniq_chain()
    print("Unique chains:", uniq_chains)

    print(ss_list)
    print(len(ss_list))
    print(len(ss_list[0]))

    SS_new_list = []

    for ss_frame in ss_list:

        ss_dict = {}
        #print(f"Secondary structure: {ss_sec}")


        for chain, seq in zip(uniq_chains, ss_frame):
            #print(f"Chain: {chain}, Sequence: {seq}")
            new_chain = ""
            for letter in chain:
                if letter != '\x00':
                    new_chain += letter
            ss_dict[new_chain] = seq
        SS_new_list.append(ss_dict)

    print(SS_new_list)
    return SS_new_list