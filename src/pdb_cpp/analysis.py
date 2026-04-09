#!/usr/bin/env python3
# coding: utf-8

import itertools
import logging

import numpy as np

from .core import align_seq_based, coor_align, get_common_atoms, rmsd as core_rmsd
from .select import remove_incomplete_backbone_residues

__all__ = ["rmsd", "interface_rmsd", "native_contact", "dockQ", "dockQ_multimer"]


logger = logging.getLogger(__name__)


def _sequence_identity(seq1, seq2):
    """Return position-wise sequence identity normalised by the longer sequence.

    Tries offsets of ±2 residues so that short N/C-terminal extensions in
    either sequence do not cause a spurious identity of zero.
    """
    s1 = seq1.replace("-", "").upper()
    s2 = seq2.replace("-", "").upper()
    if not s1 or not s2:
        return 0.0
    denom = max(len(s1), len(s2))
    best = 0.0
    for offset in range(-2, 3):  # -2, -1, 0, +1, +2
        a, b = (s1[offset:], s2) if offset >= 0 else (s1, s2[-offset:])
        matches = sum(x == y for x, y in zip(a, b))
        best = max(best, matches / denom)
    return best


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


def _check_interface_viable(coor, native_coor, rec_m, lig_m, rec_n, lig_n, back_atom, cutoff=10.0):
    """Cheap pre-check: return (True, "") when the interface is scoreable.

    Avoids the expensive ``dockQ`` call (and its logged warnings) when it
    would obviously fail because:
    - a model chain has no protein backbone atoms, or
    - the native structure has no interface residues between the two chains.
    """
    back_list = " ".join(back_atom)
    model_rec = coor.select_atoms(f"protein and chain {rec_m} and name {back_list}")
    model_lig = coor.select_atoms(f"protein and chain {lig_m} and name {back_list}")
    if len(np.unique(model_rec.models[0].uniq_resid)) == 0:
        return False, f"no backbone atoms in model chain {rec_m}"
    if len(np.unique(model_lig.models[0].uniq_resid)) == 0:
        return False, f"no backbone atoms in model chain {lig_m}"
    nat_interface = native_coor.select_atoms(
        f"chain {rec_n} and within {cutoff} of chain {lig_n}"
    )
    if len(np.unique(nat_interface.models[0].uniq_resid)) == 0:
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

    Returns
    -------
    list[float]
        Interface RMSD values for each model. Returns ``None`` entries when
        no interface residues are found.

    Notes
    -----
    Both coordinate objects must contain equivalent residues and chain order
    for the selected interface residues.
    """

    if back_atom is None:
        back_atom = ["CA", "N", "C", "O"]

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

    coor_align(coor, coor_native, index, index_native, frame_ref=0)

    return rmsd(coor, coor_native, index_list=[index, index_native])


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
        f"protein and not altloc B C D and chain {' '.join(rec_chains + lig_chains)}"
    )
    clean_native_coor = native_coor.select_atoms(
        f"protein and not altloc B C D and chain {' '.join(native_rec_chains + native_lig_chains)}"
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

    irmsd_list = interface_rmsd(
        interface_coor,
        interface_native_coor,
        native_rec_chains,
        native_lig_chains,
        cutoff=10.0,
        back_atom=back_atom,
    )
    logger.info("Interface   RMSD: %s A", irmsd_list[0])

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

    clashes_list = _count_clashes(interface_coor, rec_chains, lig_chains)
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

    if chain_map is None:
        native_chains_ordered = sorted(native_seq.keys())
        model_chains_list = list(model_seq.keys())

        # Build sequence clusters: for each native chain, find compatible model chains.
        clusters = {}
        for nat_chain in native_chains_ordered:
            nat_s = native_seq[nat_chain].replace("-", "").upper()
            candidates = [
                mc for mc in model_chains_list
                if _sequence_identity(nat_s, model_seq[mc].replace("-", "").upper()) >= 0.9
            ]
            if not candidates:
                # Fall back to the single best-matching model chain.
                candidates = [max(
                    model_chains_list,
                    key=lambda mc, ns=nat_s: _sequence_identity(
                        ns, model_seq[mc].replace("-", "").upper()
                    ),
                )]
            clusters[nat_chain] = candidates

        def _gen_maps(nat_chains, used):
            if not nat_chains:
                yield {}
                return
            nat = nat_chains[0]
            for mod in clusters[nat]:
                if mod not in used:
                    for rest in _gen_maps(nat_chains[1:], used | {mod}):
                        yield {nat: mod, **rest}

        def _score_map(candidate_map):
            try:
                trial = dockQ_multimer(
                    coor, native_coor, chain_map=candidate_map, back_atom=back_atom
                )
                valid_ifaces = [
                    v for v in trial["interfaces"].values()
                    if v is not None and any(ir is not None for ir in v["iRMS"])
                ]
                return sum(v["DockQ"][0] for v in valid_ifaces)
            except Exception as exc:
                logger.debug("Candidate chain map %s failed: %s", candidate_map, exc)
                return -1.0

        # Try the identity mapping first (native chain X → model chain X when
        # available, otherwise the first sequence-compatible candidate).  This
        # covers the common case cheaply and gives a good baseline for pruning.
        identity_map = {
            nat: (nat if nat in clusters[nat] else clusters[nat][0])
            for nat in native_chains_ordered
        }
        best_total = _score_map(identity_map)
        best_map = identity_map if best_total >= 0 else None
        logger.info("Identity mapping %s scores %.3f", identity_map, best_total)

        # Count total permutations; if only one exists the loop below is a no-op.
        _MAX_PERMUTATIONS = 720
        n_tried = 0
        for candidate_map in _gen_maps(native_chains_ordered, set()):
            if candidate_map == identity_map:
                continue  # already scored above
            n_tried += 1
            if n_tried > _MAX_PERMUTATIONS:
                logger.warning(
                    "Stopping chain-map search after %d permutations; "
                    "provide an explicit chain_map for an exhaustive search.",
                    _MAX_PERMUTATIONS,
                )
                break
            total = _score_map(candidate_map)
            if total > best_total:
                best_total = total
                best_map = candidate_map

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

        viable, reason = _check_interface_viable(
            coor, native_coor, rec_m, lig_m, rec_n, lig_n, back_atom
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

    valid = [
        v for v in interfaces.values()
        if v is not None and any(ir is not None for ir in v["iRMS"])
    ]
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
