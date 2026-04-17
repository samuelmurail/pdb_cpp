#!/usr/bin/env python3
# coding: utf-8

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geom import distance_matrix

__all__ = ["SaltBridge", "salt_bridges"]


@dataclass(frozen=True)
class SaltBridge:
    cation_resid: int
    cation_resname: str
    cation_chain: str
    cation_name: str
    cation_xyz: tuple[float, float, float]
    anion_resid: int
    anion_resname: str
    anion_chain: str
    anion_name: str
    anion_xyz: tuple[float, float, float]
    distance: float


POSITIVE_ATOMS = {
    "ARG": {"NE", "NH1", "NH2"},
    "LYS": {"NZ"},
    "HIP": {"ND1", "NE2"},
    "HSP": {"ND1", "NE2"},
}


NUCLEIC_RESNAMES = {"A", "U", "G", "C", "DA", "DT", "DC", "DG"}
NUCLEIC_ANION_ATOMS = {"OP1", "OP2", "O1P", "O2P"}


NEGATIVE_ATOMS = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
    "A": NUCLEIC_ANION_ATOMS,
    "U": NUCLEIC_ANION_ATOMS,
    "G": NUCLEIC_ANION_ATOMS,
    "C": NUCLEIC_ANION_ATOMS,
    "DA": NUCLEIC_ANION_ATOMS,
    "DT": NUCLEIC_ANION_ATOMS,
    "DC": NUCLEIC_ANION_ATOMS,
    "DG": NUCLEIC_ANION_ATOMS,
}


def _char_array_to_str_list(array_like):
    out = []
    for item in array_like:
        value = ""
        for ch in item:
            if ch != "\x00" and ch != " ":
                value += ch
        out.append(value)
    return out


def _get_model_strings(model, getter_name: str) -> list[str]:
    return _char_array_to_str_list(getattr(model, getter_name)())


def _collect_charged_atoms(model, charge_table: dict[str, set[str]]) -> list[dict]:
    names = _get_model_strings(model, "get_name")
    resnames = _get_model_strings(model, "get_resname")
    chains = _get_model_strings(model, "get_chain")
    uniq_resids = list(model.get_uniqresid())
    xs = list(model.get_x())
    ys = list(model.get_y())
    zs = list(model.get_z())

    atoms = []
    for idx, (name, resname) in enumerate(zip(names, resnames)):
        allowed = charge_table.get(resname)
        if not allowed or name not in allowed:
            continue
        atoms.append(
            {
                "resid": int(uniq_resids[idx]),
                "resname": resname,
                "chain": chains[idx],
                "name": name,
                "xyz": (float(xs[idx]), float(ys[idx]), float(zs[idx])),
            }
        )
    return atoms


def salt_bridges(
    coor,
    cation_sel: str = "protein",
    anion_sel: str = "protein",
    cutoff: float = 4.0,
):
    """Identify salt bridges between two selections for every model in *coor*.

    Salt bridges are detected between explicitly typed cationic and anionic
    heavy atoms using a simple distance cutoff.

    Parameters
    ----------
    coor : Coor
        Coordinate object (one or more models / frames).
    cation_sel : str, optional
        Atom selection string for the cationic side (default: ``"protein"``).
    anion_sel : str, optional
        Atom selection string for the anionic side (default: ``"protein"``).
    cutoff : float, optional
        Maximum cation-anion heavy-atom distance in Å (default 4.0).

    Returns
    -------
    list[list[SaltBridge]]
        One list of salt bridges per model frame.

    Notes
    -----
    Canonical nucleic acids contribute phosphate anions but no cationic groups,
    so nucleic-nucleic systems typically return no salt bridges unless modified
    positively charged residues are present.
    """
    results = []
    for frame_idx, _ in enumerate(coor.models):
        cation_model = coor.select_atoms(cation_sel, frame_idx).models[frame_idx]
        anion_model = coor.select_atoms(anion_sel, frame_idx).models[frame_idx]

        cations = _collect_charged_atoms(cation_model, POSITIVE_ATOMS)
        anions = _collect_charged_atoms(anion_model, NEGATIVE_ATOMS)

        if not cations or not anions:
            results.append([])
            continue

        cation_xyz = np.asarray([atom["xyz"] for atom in cations], dtype=np.float32)
        anion_xyz = np.asarray([atom["xyz"] for atom in anions], dtype=np.float32)
        dmat = distance_matrix(cation_xyz, anion_xyz)

        bridges = []
        seen = set()
        for cat_idx, an_idx in np.argwhere(dmat <= cutoff):
            cat = cations[int(cat_idx)]
            an = anions[int(an_idx)]
            if cat["chain"] == an["chain"] and cat["resid"] == an["resid"]:
                continue
            key = (
                cat["chain"],
                cat["resid"],
                cat["name"],
                an["chain"],
                an["resid"],
                an["name"],
            )
            if key in seen:
                continue
            seen.add(key)
            bridges.append(
                SaltBridge(
                    cation_resid=cat["resid"],
                    cation_resname=cat["resname"],
                    cation_chain=cat["chain"],
                    cation_name=cat["name"],
                    cation_xyz=cat["xyz"],
                    anion_resid=an["resid"],
                    anion_resname=an["resname"],
                    anion_chain=an["chain"],
                    anion_name=an["name"],
                    anion_xyz=an["xyz"],
                    distance=float(dmat[int(cat_idx), int(an_idx)]),
                )
            )
        results.append(bridges)
    return results