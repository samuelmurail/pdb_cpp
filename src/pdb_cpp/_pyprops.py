#!/usr/bin/env python3
# coding: utf-8

import numpy as np

from .core import Coor, Model
from . import rcsb

_coor_init = Coor.__init__
def _fetch_mmcif(pdb_id):
    """Fetch and cache the asymmetric-unit mmCIF file for a PDB ID."""
    return rcsb.download_structure(pdb_id, structure="asymmetric_unit", file_format="cif")


def _coor_init_with_pdb_id(
    self,
    coor_in=None,
    pdb_id=None,
    format="",
    rcsb_structure="asymmetric_unit",
    assembly_id=1,
    cache_dir=None,
    force_download=False,
):
    """Initialize a Coor object from a file path or PDB ID.

    Parameters
    ----------
    coor_in : str, optional
        Path to a coordinate file.
    pdb_id : str, optional
        PDB ID to fetch from RCSB.
    format : str, optional
        Force a specific file format: ``'pdb'``, ``'cif'``, ``'pqr'``, or
        ``'gro'``. Defaults to ``""`` (infer from file extension).
    rcsb_structure : str, optional
        Either ``"asymmetric_unit"`` or ``"biological_assembly"``.
    assembly_id : int, optional
        Assembly identifier when ``rcsb_structure`` is
        ``"biological_assembly"``.
    cache_dir : str, optional
        Cache directory for downloaded RCSB files.
    force_download : bool, optional
        Re-download cached RCSB files.

    Returns
    -------
    None
    """
    _coor_init(self)
    if coor_in is not None and pdb_id is not None:
        raise ValueError("Provide only one of coor_in or pdb_id")
    if pdb_id is not None:
        local_path = rcsb.download_structure(
            pdb_id,
            structure=rcsb_structure,
            file_format="cif",
            assembly_id=assembly_id,
            cache_dir=cache_dir,
            force_download=force_download,
        )
        self.read(local_path, format)
    elif coor_in is not None:
        self.read(coor_in, format)


Coor.__init__ = _coor_init_with_pdb_id


def _char_array_to_str_list(array_like):
    """Convert a list of fixed-size char arrays to strings.

    Parameters
    ----------
    array_like : iterable
        Iterable of fixed-size character arrays.

    Returns
    -------
    list[str]
        List of decoded strings with padding removed.
    """
    out = []
    for item in array_like:
        value = ""
        for ch in item:
            if ch != "\x00" and ch != " ":
                value += ch
        out.append(value)
    return out


@property
def num(self):
    """Return the atom number of the selection.

    Returns
    -------
    list[int]
        Atom indices.
    """
    return self.get_num()


@property
def uniq_resid(self):
    """Return the unique residue identifiers of the selection.

    Returns
    -------
    list[int]
        Unique residue identifiers.
    """
    return self.get_uniqresid()


@property
def name(self):
    """Return the atom names of the selection.

    Returns
    -------
    list
        Atom names.
    """
    return self.get_name()


@property
def chain(self):
    """Return the chain identifiers of the selection.

    Returns
    -------
    list
        Chain identifiers.
    """
    return self.get_chain()


@property
def resname(self):
    """Return the residue names of the selection.

    Returns
    -------
    list
        Residue names.
    """
    return self.get_resname()


@property
def resid(self):
    """Return the residue IDs of the selection.

    Returns
    -------
    list[int]
        Residue IDs.
    """
    return self.get_resid()


@property
def x(self):
    """Return the x coordinates of the selection.

    Returns
    -------
    list[float]
        X coordinates.
    """
    return self.get_x()


@property
def y(self):
    """Return the y coordinates of the selection.

    Returns
    -------
    list[float]
        Y coordinates.
    """
    return self.get_y()


@property
def z(self):
    """Return the z coordinates of the selection.

    Returns
    -------
    list[float]
        Z coordinates.
    """
    return self.get_z()


@property
def beta(self):
    """Return the B-factors of the selection.

    Returns
    -------
    list[float]
        B-factors.
    """
    return self.get_beta()


@property
def occ(self):
    """Return the occupancies of the selection.

    Returns
    -------
    list[float]
        Occupancies.
    """
    return self.get_occ()


@property
def xyz(self):
    """Return the xyz coordinates as an (N, 3) array.

    Returns
    -------
    numpy.ndarray
        Coordinate array with shape (N, 3).
    """
    return np.column_stack((self.get_x(), self.get_y(), self.get_z()))


@property
def chain_str(self):
    """Return the chain identifiers as strings.

    Returns
    -------
    list[str]
        Chain identifiers.
    """
    return _char_array_to_str_list(self.get_chain())


@property
def name_str(self):
    """Return the atom names as strings.

    Returns
    -------
    list[str]
        Atom names.
    """
    return _char_array_to_str_list(self.get_name())


@property
def resname_str(self):
    """Return the residue names as strings.

    Returns
    -------
    list[str]
        Residue names.
    """
    return _char_array_to_str_list(self.get_resname())


@property
def elem_str(self):
    """Return the element symbols as strings.

    Returns
    -------
    list[str]
        Element symbols.
    """
    return _char_array_to_str_list(self.get_elem())


@property
def alterloc_str(self):
    """Return the alternate location identifiers as strings.

    Returns
    -------
    list[str]
        Alternate location identifiers.
    """
    return _char_array_to_str_list(self.get_alterloc())


@property
def insertres_str(self):
    """Return the residue insertion codes as strings.

    Returns
    -------
    list[str]
        Residue insertion codes.
    """
    return _char_array_to_str_list(self.get_insertres())


Model.num = num
Model.uniq_resid = uniq_resid
Model.name = name
Model.chain = chain
Model.resname = resname
Model.resid = resid
Model.x = x
Model.y = y
Model.z = z
Model.beta = beta
Model.occ = occ
Model.xyz = xyz
Model.chain_str = chain_str
Model.name_str = name_str
Model.resname_str = resname_str
Model.elem_str = elem_str
Model.alterloc_str = alterloc_str
Model.insertres_str = insertres_str


def _group_atom_areas_by_residue(model, atom_areas):
    """Aggregate per-atom areas into residue records in atom order."""
    chain_str = _char_array_to_str_list(model.get_chain())
    resname_str = _char_array_to_str_list(model.get_resname())
    insertres_str = _char_array_to_str_list(model.get_insertres())
    resid = model.get_resid()
    uniq_resid = model.get_uniqresid()

    residue_lookup = {}
    residue_areas = []

    for atom_index, atom_area in enumerate(atom_areas):
        key = (
            chain_str[atom_index],
            resid[atom_index],
            insertres_str[atom_index],
            uniq_resid[atom_index],
            resname_str[atom_index],
        )
        if key not in residue_lookup:
            residue_lookup[key] = len(residue_areas)
            residue_areas.append(
                {
                    "chain": chain_str[atom_index],
                    "resid": resid[atom_index],
                    "insertres": insertres_str[atom_index],
                    "uniq_resid": uniq_resid[atom_index],
                    "resname": resname_str[atom_index],
                    "area": 0.0,
                }
            )
        residue_areas[residue_lookup[key]]["area"] += float(atom_area)

    return residue_areas


_core_model_sasa = Model.sasa


def _model_sasa_with_residue(
    self,
    probe_radius=1.4,
    n_points=960,
    include_hydrogen=False,
    by_atom=False,
    by_residue=False,
):
    """Compute SASA with optional per-atom and per-residue breakdowns."""
    need_atom_areas = by_atom or by_residue
    result = _core_model_sasa(
        self,
        probe_radius=probe_radius,
        n_points=n_points,
        include_hydrogen=include_hydrogen,
        by_atom=need_atom_areas,
    )

    if by_residue:
        result["residue_areas"] = _group_atom_areas_by_residue(self, result["atom_areas"])

    if not by_atom and "atom_areas" in result:
        del result["atom_areas"]

    return result


Model._core_sasa = _core_model_sasa
Model.sasa = _model_sasa_with_residue


@property
def length(self):  # Call it length to avoid conflict with len()
    """Return the number of atoms in the selection.

    Returns
    -------
    int
        Number of atoms.
    """
    return self.size()


Model.len = length
Coor.len = length


@property
def models(self):
    """Return all models of the selection.

    Returns
    -------
    list[Model]
        List of models.
    """
    return self.get_all_Models()


@property
def model_num(self):
    """Return the number of models in the selection.

    Returns
    -------
    int
        Number of models.
    """
    return self.model_size()


Coor.models = models
Coor.model_num = model_num


@property
def resname(self):
    """Return the residue names of the selection.

    Returns
    -------
    list
        Residue names.
    """
    return self.models[self.active_model].get_resname()


@property
def resid(self):
    """Return the residue IDs of the selection.

    Returns
    -------
    list[int]
        Residue IDs.
    """
    return self.models[self.active_model].get_resid()


@property
def uniq_resid(self):
    """Return the unique residue identifiers of the selection.

    Returns
    -------
    list[int]
        Unique residue identifiers.
    """
    return self.models[self.active_model].get_uniqresid()


@property
def chain(self):
    """Return the chain identifiers of the selection.

    Returns
    -------
    list
        Chain identifiers.
    """
    return self.models[self.active_model].get_chain()


@property
def name(self):
    """Return the atom names of the selection.

    Returns
    -------
    list
        Atom names.
    """
    return self.models[self.active_model].get_name()


@property
def num(self):
    """Return the atom numbers of the selection.

    Returns
    -------
    list[int]
        Atom indices.
    """
    return self.models[self.active_model].get_num()


@property
def x(self):
    """Return the x coordinates of the selection.

    Returns
    -------
    list[float]
        X coordinates.
    """
    return self.models[self.active_model].get_x()


@property
def y(self):
    """Return the y coordinates of the selection.

    Returns
    -------
    list[float]
        Y coordinates.
    """
    return self.models[self.active_model].get_y()


@property
def z(self):
    """Return the z coordinates of the selection.

    Returns
    -------
    list[float]
        Z coordinates.
    """
    return self.models[self.active_model].get_z()


@property
def beta(self):
    """Return the B-factors of the selection.

    Returns
    -------
    list[float]
        B-factors.
    """
    return self.models[self.active_model].get_beta()


@property
def occ(self):
    """Return the occupancies of the selection.

    Returns
    -------
    list[float]
        Occupancies.
    """
    return self.models[self.active_model].get_occ()


@property
def xyz(self):
    """Return the xyz coordinates as an (N, 3) array.

    Returns
    -------
    numpy.ndarray
        Coordinate array with shape (N, 3).
    """
    model = self.models[self.active_model]
    return np.column_stack((model.get_x(), model.get_y(), model.get_z()))


@property
def chain_str(self):
    """Return the chain identifiers as strings.

    Returns
    -------
    list[str]
        Chain identifiers.
    """
    return _char_array_to_str_list(self.models[self.active_model].get_chain())


@property
def name_str(self):
    """Return the atom names as strings.

    Returns
    -------
    list[str]
        Atom names.
    """
    return _char_array_to_str_list(self.models[self.active_model].get_name())


@property
def resname_str(self):
    """Return the residue names as strings.

    Returns
    -------
    list[str]
        Residue names.
    """
    return _char_array_to_str_list(self.models[self.active_model].get_resname())


@property
def elem_str(self):
    """Return the element symbols as strings.

    Returns
    -------
    list[str]
        Element symbols.
    """
    return _char_array_to_str_list(self.models[self.active_model].get_elem())


@property
def alterloc_str(self):
    """Return the alternate location identifiers as strings.

    Returns
    -------
    list[str]
        Alternate location identifiers.
    """
    return _char_array_to_str_list(self.models[self.active_model].get_alterloc())


@property
def insertres_str(self):
    """Return the residue insertion codes as strings.

    Returns
    -------
    list[str]
        Residue insertion codes.
    """
    return _char_array_to_str_list(self.models[self.active_model].get_insertres())


Coor.resname = resname
Coor.resid = resid
Coor.chain = chain
Coor.name = name
Coor.num = num
Coor.x = x
Coor.y = y
Coor.z = z
Coor.beta = beta
Coor.occ = occ
Coor.xyz = xyz
Coor.chain_str = chain_str
Coor.name_str = name_str
Coor.resname_str = resname_str
Coor.elem_str = elem_str
Coor.alterloc_str = alterloc_str
Coor.insertres_str = insertres_str
Coor.uniq_resid = uniq_resid
