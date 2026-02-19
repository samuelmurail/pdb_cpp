#!/usr/bin/env python3
# coding: utf-8

import logging

import numpy as np

from .core import align_seq_based, coor_align, get_common_atoms
from .select import remove_incomplete_backbone_residues

__all__ = ["rmsd", "interface_rmsd", "native_contact", "dockQ"]


logger = logging.getLogger(__name__)


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

    ref_model = coor_2.models[frame_ref]
    ref_x = np.asarray(ref_model.get_x())
    ref_y = np.asarray(ref_model.get_y())
    ref_z = np.asarray(ref_model.get_z())

    index_1 = np.asarray(index_1, dtype=int)
    index_2 = np.asarray(index_2, dtype=int)

    rmsd_list = []
    for model in coor_1.models:
        x = np.asarray(model.get_x())
        y = np.asarray(model.get_y())
        z = np.asarray(model.get_z())

        dx = x[index_1] - ref_x[index_2]
        dy = y[index_1] - ref_y[index_2]
        dz = z[index_1] - ref_z[index_2]
        rmsd_list.append(float(np.sqrt((dx * dx + dy * dy + dz * dz).mean())))

    return rmsd_list


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
    }
