#!/usr/bin/env python3
# coding: utf-8

import itertools
import logging
import math
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from .core import align_seq_based, align_index_based, coor_align, get_common_atoms, rmsd as core_rmsd
from .select import remove_incomplete_backbone_residues
from .alignment import align_seq as _align_seq

__all__ = ["rmsd", "interface_rmsd", "native_contact", "dockQ", "dockQ_multimer"]


logger = logging.getLogger(__name__)


def _count_chain_combinations(clusters):
    """Return the number of valid chain-map permutations for *clusters*.

    Mirrors ``DockQ.count_chain_combinations``.  *clusters* maps each key
    (native or model chain) to the list of compatible counterpart chains.
    The formula counts the number of ways to assign one counterpart per key
    without repeats, grouping identical candidate lists together:

    .. math::
        \\prod_{\\text{unique cluster}} P(|\\text{cluster}|, k)

    where *k* is the number of keys whose candidate list equals that cluster.
    """
    cluster_tuples = [tuple(sorted(v)) for v in clusters.values()]
    n_combos = 1
    for cluster, k in Counter(cluster_tuples).items():
        n = len(cluster)
        if n < k:
            return 0  # impossible to assign without repeats
        n_combos *= math.factorial(n) // math.factorial(n - k)
    return int(n_combos)


def _sequence_identity(seq1, seq2):
    """Return sequence identity using Smith-Waterman alignment.

    Aligns the two sequences with the C++ Smith-Waterman implementation and
    returns the fraction of identical positions within the aligned region
    (positions where neither sequence has a gap).  Normalising by the aligned
    length rather than the full sequence length means that N/C-terminal
    extensions do not artificially reduce the score.
    """
    s1 = seq1.replace("-", "").upper()
    s2 = seq2.replace("-", "").upper()
    if not s1 or not s2:
        return 0.0
    a1, a2, _ = _align_seq(s1, s2)
    aligned = [(c1, c2) for c1, c2 in zip(a1, a2) if c1 != "-" and c2 != "-"]
    if not aligned:
        return 0.0
    matches = sum(c1 == c2 for c1, c2 in aligned)
    return matches / len(aligned)


def _count_clashes(coor, rec_chains, lig_chains, clash_cutoff=2.0):
    """Count (receptor, ligand) residue pairs with any atom within *clash_cutoff* Å."""
    rec_sel = " ".join(rec_chains)
    lig_sel = " ".join(lig_chains)
    clashes_list = []
    for frame_index in range(len(coor.models)):
        close = coor.select_atoms(
            f"(chain {rec_sel} and within {clash_cutoff} of chain {lig_sel}) or "
            f"(chain {lig_sel} and within {clash_cutoff} of chain {rec_sel})",
            frame=frame_index,
        )
        rec_close = close.select_atoms(
            f"chain {rec_sel} and within {clash_cutoff} of chain {lig_sel}"
        )
        count = 0
        for residue in np.unique(rec_close.uniq_resid):
            lig_close = close.select_atoms(
                f"chain {lig_sel} and within {clash_cutoff} of residue {residue}"
            )
            count += len(np.unique(lig_close.uniq_resid))
        clashes_list.append(count)
    return clashes_list


def _check_interface_viable(coor, native_coor, rec_m, lig_m, rec_n, lig_n, back_atom, cutoff=10.0, _native_has_interface=None, _model_backbone_cache=None):
    """Cheap pre-check: return (True, "") when the interface is scoreable.

    Avoids the expensive ``dockQ`` call (and its logged warnings) when it
    would obviously fail because:
    - a model chain has no protein backbone atoms, or
    - the native structure has no interface residues between the two chains.
    """
    back_list = " ".join(back_atom)
    if _model_backbone_cache is not None:
        rec_ok = _model_backbone_cache.get(rec_m)
        lig_ok = _model_backbone_cache.get(lig_m)
    else:
        rec_ok = lig_ok = None
    if rec_ok is None:
        model_rec = coor.select_atoms(f"protein and chain {rec_m} and name {back_list}")
        rec_ok = len(np.unique(model_rec.models[0].uniq_resid)) > 0
    if not rec_ok:
        return False, f"no backbone atoms in model chain {rec_m}"
    if lig_ok is None:
        model_lig = coor.select_atoms(f"protein and chain {lig_m} and name {back_list}")
        lig_ok = len(np.unique(model_lig.models[0].uniq_resid)) > 0
    if not lig_ok:
        return False, f"no backbone atoms in model chain {lig_m}"
    if _native_has_interface is None:
        nat_interface = native_coor.select_atoms(
            f"chain {rec_n} and within {cutoff} of chain {lig_n}"
        )
        _native_has_interface = len(np.unique(nat_interface.models[0].uniq_resid)) > 0
    if not _native_has_interface:
        return False, f"no interface residues between native chains {rec_n} and {lig_n}"
    return True, ""


def rmsd(coor_1, coor_2, selection="name CA", index_list=None, frame_ref=0):
    """Compute RMSD between two sets of coordinates.

    Parameters
    ----------
    coor_1 : Coor
        First set of coordinates.
    coor_2 : Coor
        Second set of coordinates.
    selection : str, optional
        Selection string used when index_list is not provided.
    index_list : list, optional
        Pair of index lists [index_1, index_2].
    frame_ref : int, optional
        Reference frame index in coor_2.

    Returns
    -------
    list[float]
        RMSD values for each model in coor_1.
    """

    if frame_ref < 0 or frame_ref >= coor_2.model_num:
        raise ValueError(
            "Reference frame index is larger than the number of frames in the reference structure"
        )

    if index_list is None:
        index_1 = coor_1.get_index_select(selection)
        index_2 = coor_2.get_index_select(selection)
    else:
        index_1 = index_list[0]
        index_2 = index_list[1]

    if len(index_1) == 0 or len(index_2) == 0:
        raise ValueError("No atoms selected for RMSD calculation")
    if len(index_1) != len(index_2):
        raise ValueError("Index lists must have the same length")

    return core_rmsd(coor_1, coor_2, list(index_1), list(index_2), frame_ref=frame_ref)


def interface_rmsd(
    coor,
    coor_native,
    rec_chains_native,
    lig_chains_native,
    cutoff=10.0,
    back_atom=None,
    index_pair=None,
):
    """Compute the interface RMSD between two models.

    The interface is defined as atoms within ``cutoff`` Angstrom of the
    opposite partner chain(s) in the native structure. The RMSD is computed
    on the selected backbone atoms after aligning the model to the native
    structure using the interface atoms.

    Parameters
    ----------
    coor : Coor
        Model coordinates.
    coor_native : Coor
        Native coordinates.
    rec_chains_native : list[str]
        Native receptor chains.
    lig_chains_native : list[str]
        Native ligand chains.
    cutoff : float, default=10.0
        Interface distance cutoff in Angstrom.
    back_atom : list[str], optional
        Backbone atom names used for RMSD (default: ``["CA", "N", "C", "O"]``).
    index_pair : tuple[list[int], list[int]], optional
        Pre-computed ``(model_indices, native_indices)`` of the interface
        backbone atoms.  When supplied the function skips the chain-selection
        and residue-mapping steps entirely, calling :func:`align_index_based`
        directly.  This correctly handles non-sequential or non-contiguous
        index lists (e.g. when model chains are numbered from 1 on every chain
        and thus have overlapping ``resid`` values).

    Returns
    -------
    list[float]
        Interface RMSD values for each model. Returns ``None`` entries when
        no interface residues are found.
    """

    if back_atom is None:
        back_atom = ["CA", "N", "C", "O"]

    # Fast path: caller provides pre-matched atom indices.
    if index_pair is not None:
        iface_model_idx, iface_native_idx = index_pair
        if not iface_model_idx:
            return [None] * len(coor.models)
        rmsds, _, _ = align_index_based(coor, coor_native, iface_model_idx, iface_native_idx)
        return rmsds

    lig_interface = coor_native.select_atoms(
        f"chain {' '.join(lig_chains_native)} and within {cutoff} of chain {' '.join(rec_chains_native)}"
    )
    rec_interface = coor_native.select_atoms(
        f"chain {' '.join(rec_chains_native)} and within {cutoff} of chain {' '.join(lig_chains_native)}"
    )

    lig_interface_residues = np.unique(lig_interface.models[0].uniq_resid)
    rec_interface_residues = np.unique(rec_interface.models[0].uniq_resid)

    interface_residues = np.concatenate(
        (lig_interface_residues, rec_interface_residues)
    )

    if len(interface_residues) == 0:
        logger.info("No interface residues found")
        return [None] * len(coor.models)

    native_unique = list(dict.fromkeys(coor_native.models[0].uniq_resid))
    coor_unique = list(dict.fromkeys(coor.models[0].uniq_resid))
    if len(native_unique) != len(coor_unique):
        raise ValueError("Interface selections do not contain the same residue count")

    native_to_coor = {nat: coor_unique[i] for i, nat in enumerate(native_unique)}
    mapped_interface = [native_to_coor[nat] for nat in interface_residues]

    native_list = " ".join(str(i) for i in interface_residues)
    coor_list = " ".join(str(i) for i in mapped_interface)
    back_list = " ".join(back_atom)

    index = coor.get_index_select(f"residue {coor_list} and name {back_list}")
    index_native = coor_native.get_index_select(
        f"residue {native_list} and name {back_list}"
    )

    if len(index) != len(index_native):
        raise ValueError(
            "The number of atoms in the interface is not the same in the two models"
        )

    rmsds, _, _ = align_index_based(coor, coor_native, index, index_native)
    return rmsds


def native_contact(
    coor,
    native_coor,
    rec_chains,
    lig_chains,
    native_rec_chains,
    native_lig_chains,
    cutoff=5.0,
    residue_id_map=None,
    native_residue_id_map=None,
):
    """Compute native and non-native contact fractions between model and native.

    The function builds the set of native receptor-ligand residue contacts
    within ``cutoff`` Angstrom in the native structure, then counts which of
    those contacts are present in the model (Fnat) and which model contacts are
    not native (Fnonnat).

    Parameters
    ----------
    coor : Coor
        Model coordinates.
    native_coor : Coor
        Native coordinates.
    rec_chains : list[str]
        Model receptor chains.
    lig_chains : list[str]
        Model ligand chains.
    native_rec_chains : list[str]
        Native receptor chains.
    native_lig_chains : list[str]
        Native ligand chains.
    cutoff : float, default=5.0
        Contact distance cutoff in Angstrom.
    residue_id_map : dict[int, int], optional
        Mapping from model residue IDs to a shared residue ID space.
    native_residue_id_map : dict[int, int], optional
        Mapping from native residue IDs to the same shared residue ID space.

    Returns
    -------
    tuple[list[float], list[float]]
        ``(fnat_list, fnonnat_list)`` for each model in ``coor``.
    """

    if residue_id_map is None:
        residue_id_map = {}
    if native_residue_id_map is None:
        native_residue_id_map = {}

    native_rec_lig_interface = native_coor.select_atoms(
        f"(chain {' '.join(native_rec_chains)} and within {cutoff} of chain {' '.join(native_lig_chains)}) or"
        f"(chain {' '.join(native_lig_chains)} and within {cutoff} of chain {' '.join(native_rec_chains)})"
    )
    native_rec_interface = native_rec_lig_interface.select_atoms(
        f"chain {' '.join(native_rec_chains)} and within {cutoff} of chain {' '.join(native_lig_chains)}"
    )

    native_contact_list = []
    native_rec_residue = native_rec_interface.models[0].uniq_resid

    for residue in np.unique(native_rec_residue):
        native_residue_common = native_residue_id_map.get(residue, residue)
        res_lig = native_rec_lig_interface.select_atoms(
            f"chain {' '.join(native_lig_chains)} and within {cutoff} of residue {residue}"
        )
        native_contact_list += [
            [native_residue_common, native_residue_id_map.get(lig_residue, lig_residue)]
            for lig_residue in np.unique(res_lig.models[0].uniq_resid)
        ]

    fnat_list = []
    fnonnat_list = []
    for frame_index, _ in enumerate(coor.models):
        rec_lig_interface = coor.select_atoms(
            f"(chain {' '.join(rec_chains)} and within {cutoff} of chain {' '.join(lig_chains)}) or "
            f"(chain {' '.join(lig_chains)} and within {cutoff} of chain {' '.join(rec_chains)})",
            frame=frame_index,
        )
        rec_interface = rec_lig_interface.select_atoms(
            f"(chain {' '.join(rec_chains)} and within {cutoff} of chain {' '.join(lig_chains)})"
        )

        model_rec_residue = np.unique(rec_interface.uniq_resid)

        native_contact_num = 0
        non_native_contact_num = 0
        model_contact_list = []

        for residue in model_rec_residue:
            model_residue_common = residue_id_map.get(residue, residue)
            residue_lig = rec_lig_interface.select_atoms(
                f"chain {' '.join(lig_chains)} and within {cutoff} of residue {residue}"
            )

            for lig_residue in np.unique(residue_lig.uniq_resid):
                contact_pair = [
                    model_residue_common,
                    residue_id_map.get(lig_residue, lig_residue),
                ]
                if contact_pair in native_contact_list:
                    native_contact_num += 1
                else:
                    non_native_contact_num += 1

                model_contact_list.append(contact_pair)

        if native_contact_num > 0:
            fnat = native_contact_num / len(native_contact_list)
        else:
            fnat = 0.0
        logger.info(
            "Fnat %.3f %d correct of %d native contacts",
            fnat,
            native_contact_num,
            len(native_contact_list),
        )

        if non_native_contact_num > 0:
            fnonnat = non_native_contact_num / len(model_contact_list)
        else:
            fnonnat = 0.0
        logger.info(
            "Fnonnat %.3f %d non-native of %d model contacts",
            fnonnat,
            non_native_contact_num,
            len(model_contact_list),
        )

        fnat_list.append(fnat)
        fnonnat_list.append(fnonnat)

    return fnat_list, fnonnat_list


def dockQ(
    coor,
    native_coor,
    rec_chains=None,
    lig_chains=None,
    native_rec_chains=None,
    native_lig_chains=None,
    back_atom=None,
    _search_mode=False,
):
    """Compute DockQ scores between a model and a native structure.

    DockQ combines interface contacts (Fnat), ligand RMSD (LRMS), and
    interface RMSD (iRMS) into a single docking quality metric. Chain roles
    are inferred by selecting the shortest chain as ligand when not provided.

    Parameters
    ----------
    coor : Coor
        Model coordinates.
    native_coor : Coor
        Native coordinates.
    rec_chains : list[str], optional
        Model receptor chains. If ``None``, uses all chains except the
        shortest chain (ligand).
    lig_chains : list[str], optional
        Model ligand chains. If ``None``, uses the shortest chain.
    native_rec_chains : list[str], optional
        Native receptor chains. If ``None``, uses all chains except the
        shortest chain (ligand).
    native_lig_chains : list[str], optional
        Native ligand chains. If ``None``, uses the shortest chain.
    back_atom : list[str], optional
        Backbone atom names used for alignment and RMSD calculations.

    Returns
    -------
    dict
        Dictionary with keys ``Fnat``, ``Fnonnat``, ``rRMS``, ``iRMS``,
        ``LRMS``, and ``DockQ``, each containing lists per model.

    Notes
    -----
    This implementation mirrors the pdb_numpy DockQ pipeline and relies on
    sequence-based alignment for receptor superposition before computing the
    ligand RMSD and interface metrics.
    """

    if back_atom is None:
        back_atom = ["CA", "N", "C", "O"]

    model_seq = coor.get_aa_seq()
    native_seq = native_coor.get_aa_seq()

    if lig_chains is None:
        lig_chains = [
            min(model_seq.items(), key=lambda x: len(x[1].replace("-", "")))[0]
        ]
    logger.info("Model ligand chains : %s", " ".join(lig_chains))
    if rec_chains is None:
        rec_chains = [chain for chain in model_seq if chain not in lig_chains]
    logger.info("Model receptor chains : %s", " ".join(rec_chains))

    if native_lig_chains is None:
        native_lig_chains = [
            min(native_seq.items(), key=lambda x: len(x[1].replace("-", "")))[0]
        ]
    logger.info("Native ligand chains : %s", " ".join(native_lig_chains))
    if native_rec_chains is None:
        native_rec_chains = [
            chain for chain in native_seq if chain not in native_lig_chains
        ]
    logger.info("Native receptor chains : %s", " ".join(native_rec_chains))

    clean_coor = coor.select_atoms(
        f"protein and not altloc B C D and not name H* and chain {' '.join(rec_chains + lig_chains)}"
    )
    clean_native_coor = native_coor.select_atoms(
        f"protein and not altloc B C D and not name H* and chain {' '.join(native_rec_chains + native_lig_chains)}"
    )

    clean_coor = remove_incomplete_backbone_residues(clean_coor, back_atom)
    clean_native_coor = remove_incomplete_backbone_residues(
        clean_native_coor, back_atom
    )

    rmsd_prot_list, align_rec_index, align_rec_native_index = align_seq_based(
        clean_coor,
        clean_native_coor,
        chain_1=rec_chains,
        chain_2=native_rec_chains,
        back_names=back_atom,
    )
    logger.info("Receptor RMSD: %.3f A", rmsd_prot_list[0])
    logger.info(
        "Found %d residues in common (receptor)",
        len(align_rec_index) // len(back_atom),
    )

    lig_index, lig_native_index = get_common_atoms(
        clean_coor,
        clean_native_coor,
        chain_1=lig_chains,
        chain_2=native_lig_chains,
        back_names=back_atom,
        matrix_file="",
    )
    lrmsd_list = rmsd(
        clean_coor, clean_native_coor, index_list=[lig_index, lig_native_index]
    )
    logger.info("Ligand   RMSD: %.3f A", lrmsd_list[0])
    logger.info(
        "Found %d residues in common (ligand)",
        len(lig_index) // len(back_atom),
    )

    # Fast path for chain-map search: only LRMS matters for ranking.
    # Skips native_contact (expensive Python loops), interface_rmsd, and all
    # intermediate select_atoms / within queries.  The returned pseudo-DockQ is
    # monotonically related to real DockQ for correct vs. incorrect mappings.
    if _search_mode:
        def _scale(rms, d):
            return 0.0 if rms is None else 1.0 / (1.0 + (rms / d) ** 2)
        pseudo_dq = [_scale(lrms, 8.5) / 3.0 for lrms in lrmsd_list]
        return {
            "Fnat": [0.0] * len(lrmsd_list),
            "Fnonnat": [0.0] * len(lrmsd_list),
            "rRMS": rmsd_prot_list,
            "iRMS": [None] * len(lrmsd_list),
            "LRMS": lrmsd_list,
            "DockQ": pseudo_dq,
            "clashes": [0] * len(lrmsd_list),
        }

    coor_residue = np.asarray(clean_coor.models[0].uniq_resid)[
        np.asarray(align_rec_index + lig_index, dtype=int)
    ]
    native_residue = np.asarray(clean_native_coor.models[0].uniq_resid)[
        np.asarray(align_rec_native_index + lig_native_index, dtype=int)
    ]

    coor_residue_unique = np.unique(coor_residue)
    native_residue_unique = np.unique(native_residue)

    if len(coor_residue_unique) != len(native_residue_unique):
        raise ValueError("Residue alignment mismatch between model and native")

    interface_coor = clean_coor.select_atoms(
        f"residue {' '.join(str(i) for i in coor_residue_unique)}"
    )
    interface_native_coor = clean_native_coor.select_atoms(
        f"residue {' '.join(str(i) for i in native_residue_unique)}"
    )

    residue_id_map = {}
    native_residue_id_map = {}
    common_residue_ids = {}
    common_id = 1
    for model_residue_id, native_residue_id in zip(coor_residue, native_residue):
        if native_residue_id not in common_residue_ids:
            common_residue_ids[native_residue_id] = common_id
            common_id += 1

        mapped_id = common_residue_ids[native_residue_id]
        residue_id_map[model_residue_id] = mapped_id
        native_residue_id_map[native_residue_id] = mapped_id

    # Build interface atom index pairs from the sequence-aligned indices.
    # Filtering to interface residues here (native side within 10 Å of partner)
    # lets us pass them directly to interface_rmsd via index_pair, which calls
    # align_index_based.  This is correct even when model chains restart
    # residue numbering from 1 (e.g. CASP models) because we work with atom
    # indices, not residue selections.
    nat_lig_iface = clean_native_coor.select_atoms(
        f"chain {' '.join(native_lig_chains)} "
        f"and within 10.0 of chain {' '.join(native_rec_chains)}"
    )
    nat_rec_iface = clean_native_coor.select_atoms(
        f"chain {' '.join(native_rec_chains)} "
        f"and within 10.0 of chain {' '.join(native_lig_chains)}"
    )
    iface_native_uids = set(
        list(nat_lig_iface.models[0].uniq_resid)
        + list(nat_rec_iface.models[0].uniq_resid)
    )

    native_all_uids = np.asarray(clean_native_coor.models[0].uniq_resid)
    all_model_idx = align_rec_index + lig_index
    all_native_idx = align_rec_native_index + lig_native_index
    iface_model_idx = [
        mi for mi, ni in zip(all_model_idx, all_native_idx)
        if native_all_uids[ni] in iface_native_uids
    ]
    iface_native_idx = [
        ni for mi, ni in zip(all_model_idx, all_native_idx)
        if native_all_uids[ni] in iface_native_uids
    ]

    # Fnat is computed on clean_coor (before interface_rmsd aligns it in-place).
    fnat_list, fnonnat_list = native_contact(
        interface_coor,
        interface_native_coor,
        rec_chains,
        lig_chains,
        native_rec_chains,
        native_lig_chains,
        cutoff=5.0,
        residue_id_map=residue_id_map,
        native_residue_id_map=native_residue_id_map,
    )
    logger.info("Fnat: %.3f      Fnonnat: %.3f", fnat_list[0], fnonnat_list[0])

    if _search_mode:
        clashes_list = [0] * len(fnat_list)
    else:
        clashes_list = _count_clashes(interface_coor, rec_chains, lig_chains)

    irmsd_list = interface_rmsd(
        clean_coor,
        clean_native_coor,
        native_rec_chains,
        native_lig_chains,
        cutoff=10.0,
        back_atom=back_atom,
        index_pair=(iface_model_idx, iface_native_idx) if iface_model_idx else None,
    )
    logger.info("Interface   RMSD: %s A", irmsd_list[0])
    logger.info("Clashes: %d", clashes_list[0])

    def scale_rms(rms, d):
        if rms is None:
            return 0.0
        return 1.0 / (1 + (rms / d) ** 2)

    d1 = 8.5
    d2 = 1.5

    dockq_list = [
        (fnat + scale_rms(lrmsd, d1) + scale_rms(irmsd, d2)) / 3
        for fnat, lrmsd, irmsd in zip(fnat_list, lrmsd_list, irmsd_list)
    ]

    logger.info("DockQ Score pdb_cpp: %.3f", dockq_list[0])

    return {
        "Fnat": fnat_list,
        "Fnonnat": fnonnat_list,
        "rRMS": rmsd_prot_list,
        "iRMS": irmsd_list,
        "LRMS": lrmsd_list,
        "DockQ": dockq_list,
        "clashes": clashes_list,
    }


def dockQ_multimer(
    coor,
    native_coor,
    chain_map=None,
    back_atom=None,
    n_cpu=1,
    _search_mode=False,
    _native_iface_cache=None,
    _model_backbone_cache=None,
):
    """Compute DockQ over all pairwise native chain interfaces (multimer).

    Scores every :math:`\\binom{n}{2}` interface between the *n* native chains
    and returns per-interface DockQ metrics as well as *GlobalDockQ* (the
    average DockQ over all interfaces), mirroring the DockQ v2 multimer output.

    Parameters
    ----------
    coor : Coor
        Model coordinates.
    native_coor : Coor
        Native coordinates.
    chain_map : dict[str, str], optional
        Mapping from **native** chain IDs to **model** chain IDs.  If ``None``
        the function assumes that both structures share the same chain names
        and builds an identity mapping for chains present in both.
    back_atom : list[str], optional
        Backbone atom names used for alignment and RMSD calculations.

    Returns
    -------
    dict
        Dictionary with two keys:

        ``"interfaces"``
            ``dict[(native_ch1, native_ch2), result]`` where *result* is the
            dict returned by :func:`dockQ` for that pair (or ``None`` when the
            interface could not be scored).
        ``"GlobalDockQ"``
            ``list[float]`` — average DockQ over all valid interfaces,
            one value per model frame.

    Notes
    -----
    The larger native chain of each pair is used as receptor to match the
    DockQ v2 convention.  For chains of equal length the pair is presented in
    the order they appear in the ``chain_map`` iteration order.
    """
    if back_atom is None:
        back_atom = ["CA", "N", "C", "O"]

    native_seq = native_coor.get_aa_seq()
    model_seq = coor.get_aa_seq()

    # These are populated inside the auto-mapping block and reused in the
    # final per-interface scoring pass below, avoiding redundant select_atoms.
    _native_iface_cache_local = _native_iface_cache
    _model_backbone_cache_local = _model_backbone_cache

    if chain_map is None:
        native_chains_ordered = sorted(native_seq.keys())
        model_chains_list = sorted(model_seq.keys())

        _MAX_PERMUTATIONS = 720
        _all_native = sorted(native_seq.keys())

        # ── Step 1: native interface presence (per pair, independent of mapping) ──
        _native_iface_cache_local = {}
        for _n1, _n2 in itertools.combinations(_all_native, 2):
            _sel = native_coor.select_atoms(f"chain {_n1} and within 10.0 of chain {_n2}")
            _native_iface_cache_local[(_n1, _n2)] = (
                len(np.unique(_sel.models[0].uniq_resid)) > 0
            )

        # ── Step 2: per model-chain backbone presence ──────────────────────────
        _back_list = " ".join(back_atom)
        _model_backbone_cache_local = {}
        for _mc in model_chains_list:
            _sel = coor.select_atoms(f"protein and chain {_mc} and name {_back_list}")
            _model_backbone_cache_local[_mc] = len(np.unique(_sel.models[0].uniq_resid)) > 0

        # ── Step 3: build sequence clusters ───────────────────────────────────
        if len(model_chains_list) < len(native_chains_ordered):
            # Inverse case: clusters maps model → [compatible native chains]
            clusters = {}
            for mod_chain in model_chains_list:
                mod_s = model_seq[mod_chain].replace("-", "").upper()
                candidates = [
                    nc for nc in native_chains_ordered
                    if _sequence_identity(mod_s, native_seq[nc].replace("-", "").upper()) >= 0.9
                ]
                if not candidates:
                    candidates = [max(
                        native_chains_ordered,
                        key=lambda nc, ms=mod_s: _sequence_identity(
                            ms, native_seq[nc].replace("-", "").upper()
                        ),
                    )]
                clusters[mod_chain] = candidates
            # Reverse lookup: native chain → model chains that could be assigned to it
            _model_cands_for_nat = {}
            for _mc, _nats in clusters.items():
                for _nc in _nats:
                    _model_cands_for_nat.setdefault(_nc, []).append(_mc)
        else:
            # Standard case: clusters maps native → [compatible model chains]
            clusters = {}
            for nat_chain in native_chains_ordered:
                nat_s = native_seq[nat_chain].replace("-", "").upper()
                candidates = [
                    mc for mc in model_chains_list
                    if _sequence_identity(nat_s, model_seq[mc].replace("-", "").upper()) >= 0.9
                ]
                if not candidates:
                    candidates = [max(
                        model_chains_list,
                        key=lambda mc, ns=nat_s: _sequence_identity(
                            ns, model_seq[mc].replace("-", "").upper()
                        ),
                    )]
                clusters[nat_chain] = candidates
            _model_cands_for_nat = clusters  # already nat → [model chains]

        # ── Step 4: precompute pairwise dockQ scores (search mode) ────────────
        # For each native interface pair (rec_n, lig_n) and every candidate
        # (rec_m, lig_m) that a permutation could assign to it, run dockQ once
        # in _search_mode (pure C++ LRMS path).  Scoring a permutation then
        # costs only O(n_interfaces) dict lookups instead of O(n_interfaces)
        # full dockQ calls, reducing 720 × 21 calls to ≲ C(n_model,2) calls.
        _native_pairs_with_iface = [
            (_n1, _n2)
            for _n1, _n2 in itertools.combinations(_all_native, 2)
            if _native_iface_cache_local.get((_n1, _n2), False)
        ]
        _pairwise_cache: dict = {}
        for _n1, _n2 in _native_pairs_with_iface:
            _n1_len = len(native_seq.get(_n1, "").replace("-", ""))
            _n2_len = len(native_seq.get(_n2, "").replace("-", ""))
            _rec_n, _lig_n = (_n1, _n2) if _n1_len >= _n2_len else (_n2, _n1)
            for _rec_m in _model_cands_for_nat.get(_rec_n, []):
                if not _model_backbone_cache_local.get(_rec_m, True):
                    continue
                for _lig_m in _model_cands_for_nat.get(_lig_n, []):
                    if _rec_m == _lig_m:
                        continue
                    if not _model_backbone_cache_local.get(_lig_m, True):
                        continue
                    _key = (_rec_m, _lig_m, _rec_n, _lig_n)
                    if _key in _pairwise_cache:
                        continue
                    try:
                        _res = dockQ(
                            coor, native_coor,
                            rec_chains=[_rec_m], lig_chains=[_lig_m],
                            native_rec_chains=[_rec_n], native_lig_chains=[_lig_n],
                            back_atom=back_atom, _search_mode=True,
                        )
                        _pairwise_cache[_key] = _res["DockQ"][0]
                    except Exception as _exc:
                        logger.debug(
                            "Pairwise precomp (%s,%s)→(%s,%s) failed: %s",
                            _rec_m, _lig_m, _rec_n, _lig_n, _exc,
                        )
                        _pairwise_cache[_key] = 0.0
        logger.info(
            "Precomputed %d pairwise scores for %d native interface pairs",
            len(_pairwise_cache), len(_native_pairs_with_iface),
        )

        # ── Step 5: score a candidate map using the precomputed matrix ─────────
        def _score_map(candidate_map):
            total = 0.0
            for _n1, _n2 in _native_pairs_with_iface:
                _m1 = candidate_map.get(_n1)
                _m2 = candidate_map.get(_n2)
                if _m1 is None or _m2 is None or _m1 == _m2:
                    continue
                _n1_len = len(native_seq.get(_n1, "").replace("-", ""))
                _n2_len = len(native_seq.get(_n2, "").replace("-", ""))
                if _n1_len >= _n2_len:
                    _rec_m, _lig_m, _rec_n, _lig_n = _m1, _m2, _n1, _n2
                else:
                    _rec_m, _lig_m, _rec_n, _lig_n = _m2, _m1, _n2, _n1
                total += _pairwise_cache.get((_rec_m, _lig_m, _rec_n, _lig_n), 0.0)
            return total

        def _find_best_map(all_maps, initial_map=None):
            best_total = -1.0
            best_map = None
            if initial_map is not None:
                best_total = _score_map(initial_map)
                best_map = initial_map if best_total >= 0 else None
                logger.info("Identity mapping %s scores %.3f", initial_map, best_total)
            maps_list = [m for m in all_maps if m != initial_map]
            if len(maps_list) > _MAX_PERMUTATIONS:
                logger.warning(
                    "Stopping chain-map search after %d permutations (found %d); "
                    "provide an explicit chain_map for an exhaustive search.",
                    _MAX_PERMUTATIONS, len(maps_list),
                )
                maps_list = maps_list[:_MAX_PERMUTATIONS]
            if not maps_list:
                return best_map, best_total
            for cmap in maps_list:
                total = _score_map(cmap)
                if total > best_total:
                    best_total = total
                    best_map = cmap
            return best_map, best_total

        # ── Step 6: generate candidate maps and pick the best ─────────────────
        if len(model_chains_list) < len(native_chains_ordered):
            def _gen_maps_inv(mod_chains, used):
                """Yield model→native dicts; each native chain used at most once."""
                if not mod_chains:
                    yield {}
                    return
                mod = mod_chains[0]
                for nat in clusters[mod]:
                    if nat not in used:
                        for rest in _gen_maps_inv(mod_chains[1:], used | {nat}):
                            yield {mod: nat, **rest}

            n_combos = _count_chain_combinations(clusters)
            logger.info(
                "Chain-map search (inverse): %d permutation(s) to evaluate", n_combos
            )
            all_n2m = [
                {v: k for k, v in m2n.items()}
                for m2n in _gen_maps_inv(model_chains_list, set())
            ]
            best_map, best_total = _find_best_map(all_n2m)

        else:
            def _gen_maps(nat_chains, used):
                if not nat_chains:
                    yield {}
                    return
                nat = nat_chains[0]
                for mod in clusters[nat]:
                    if mod not in used:
                        for rest in _gen_maps(nat_chains[1:], used | {mod}):
                            yield {nat: mod, **rest}

            n_combos = _count_chain_combinations(clusters)
            logger.info(
                "Chain-map search: %d permutation(s) to evaluate", n_combos
            )
            identity_map = {
                nat: (nat if nat in clusters[nat] else clusters[nat][0])
                for nat in native_chains_ordered
            }
            all_nat_maps = list(_gen_maps(native_chains_ordered, set()))
            best_map, best_total = _find_best_map(all_nat_maps, initial_map=identity_map)

        if best_map is None:
            best_map = {ch: ch for ch in native_seq if ch in model_seq}
            logger.warning("No valid chain mapping found; falling back to identity mapping.")
        chain_map = best_map
        logger.info("Optimal chain mapping: %s (sum DockQ: %.3f)", chain_map, best_total)

    native_chains = sorted(chain_map.keys())
    interfaces = {}

    for n1, n2 in itertools.combinations(native_chains, 2):
        m1 = chain_map[n1]
        m2 = chain_map[n2]

        # Skip if both native chains map to the same model chain (can happen
        # when model has fewer chains than native, e.g. a subcomplex model vs.
        # a full biological assembly). There is no inter-chain interface to score.
        if m1 == m2:
            logger.info(
                "Skipping interface (%s, %s): both map to the same model chain %s",
                n1, n2, m1,
            )
            interfaces[(n1, n2)] = None
            continue

        # Assign the larger native chain as receptor (DockQ v2 convention)
        n1_len = len(native_seq.get(n1, "").replace("-", ""))
        n2_len = len(native_seq.get(n2, "").replace("-", ""))
        if n1_len >= n2_len:
            rec_n, lig_n, rec_m, lig_m = n1, n2, m1, m2
        else:
            rec_n, lig_n, rec_m, lig_m = n2, n1, m2, m1

        logger.info(
            "Scoring interface native(%s, %s) -> model(%s, %s)",
            rec_n, lig_n, rec_m, lig_m,
        )

        # Use precomputed caches — prefer locally computed ones (set above)
        # which are available when we auto-detected the chain map.
        _eff_nat_cache = _native_iface_cache_local
        _cached_nat = (
            _eff_nat_cache.get((n1, n2), _eff_nat_cache.get((n2, n1)))
            if _eff_nat_cache is not None else None
        )
        viable, reason = _check_interface_viable(
            coor, native_coor, rec_m, lig_m, rec_n, lig_n, back_atom,
            _native_has_interface=_cached_nat,
            _model_backbone_cache=_model_backbone_cache_local,
        )
        if not viable:
            logger.info("Skipping interface (%s, %s): %s", n1, n2, reason)
            interfaces[(n1, n2)] = None
            continue

        try:
            result = dockQ(
                coor,
                native_coor,
                rec_chains=[rec_m],
                lig_chains=[lig_m],
                native_rec_chains=[rec_n],
                native_lig_chains=[lig_n],
                back_atom=back_atom,
                _search_mode=_search_mode,
            )
            interfaces[(n1, n2)] = {
                **result,
                "model_rec_chain": rec_m,
                "model_lig_chain": lig_m,
            }
        except Exception as exc:
            logger.info(
                "Could not compute DockQ for interface (%s, %s): %s", n1, n2, exc
            )
            interfaces[(n1, n2)] = None

    # iRMS=None just causes scale_rms to return 0.0 in the DockQ formula;
    # the score is still well-defined, so we only drop interfaces where dockQ
    # raised an exception (v is None).
    valid = [v for v in interfaces.values() if v is not None]
    if valid:
        n_frames = len(valid[0]["DockQ"])
        global_dockq = [
            sum(v["DockQ"][i] for v in valid) / len(valid)
            for i in range(n_frames)
        ]
    else:
        logger.warning("No valid interfaces found; GlobalDockQ is undefined.")
        global_dockq = []

    logger.info("GlobalDockQ (avg over %d interfaces): %.3f", len(valid), global_dockq[0] if global_dockq else float("nan"))

    return {
        "interfaces": interfaces,
        "GlobalDockQ": global_dockq,
        "chain_map": chain_map,
    }
