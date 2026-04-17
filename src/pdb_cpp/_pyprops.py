#!/usr/bin/env python3
# coding: utf-8

import numpy as np

from .core import Coor, Model
from . import rcsb


class _FieldProxy:
    """Proxy for a per-atom field that supports indexed get/set.

    Holds a reference to the original C++ object (Coor or Model) and
    calls the named getter/setter methods directly, so writes always
    propagate back into the underlying data.

    Supports::

        obj.beta[indices] = 42.0      # broadcast scalar
        obj.x[indices] = array        # assign array
        obj.chain[indices] = "B"      # broadcast string
    """

    __slots__ = ("_target", "_getter", "_setter")

    def __init__(self, target, getter: str, setter: str):
        object.__setattr__(self, "_target", target)
        object.__setattr__(self, "_getter", getter)
        object.__setattr__(self, "_setter", setter)

    def _raw(self):
        return getattr(self._target, self._getter)()

    def _get_array(self):
        return np.asarray(self._raw(), dtype=float)

    def __getitem__(self, index):
        raw = self._raw()
        if isinstance(raw, (list, tuple)):
            if np.isscalar(index):
                return raw[index]
            return [raw[i] for i in index]
        return np.asarray(raw)[index]

    def __setitem__(self, index, value):
        setter = getattr(self._target, self._setter)
        indices = [int(index)] if np.isscalar(index) else [int(i) for i in index]
        if isinstance(value, str) or (np.isscalar(value) and not isinstance(value, (list, tuple))):
            for i in indices:
                setter(i, value)
        else:
            for i, v in zip(indices, value):
                setter(i, v)

    def __len__(self):
        return len(self._raw())

    def __iter__(self):
        return iter(self._raw())

    def __repr__(self):
        return repr(list(self._raw()))

    def __eq__(self, other):
        left = list(self._raw())
        if isinstance(other, _FieldProxy):
            return left == list(other._raw())
        return left == other

    def __ne__(self, other):
        return not self.__eq__(other)

    def __array__(self, dtype=None):
        arr = self._get_array()
        return arr.astype(dtype) if dtype is not None else arr


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
    return _FieldProxy(self, 'get_num', 'set_num')


@property
def uniq_resid(self):
    return _FieldProxy(self, 'get_uniqresid', 'set_uniqresid')


@property
def name(self):
    return _FieldProxy(self, 'get_name', 'set_name')


@property
def chain(self):
    return _FieldProxy(self, 'get_chain', 'set_chain')


@property
def resname(self):
    return _FieldProxy(self, 'get_resname', 'set_resname')


@property
def resid(self):
    return _FieldProxy(self, 'get_resid', 'set_resid')


@property
def x(self):
    return _FieldProxy(self, 'get_x', 'set_x')


@property
def y(self):
    return _FieldProxy(self, 'get_y', 'set_y')


@property
def z(self):
    return _FieldProxy(self, 'get_z', 'set_z')


@property
def beta(self):
    return _FieldProxy(self, 'get_beta', 'set_beta')


@property
def occ(self):
    return _FieldProxy(self, 'get_occ', 'set_occ')


@property
def xyz(self):
    return np.column_stack((self.get_x(), self.get_y(), self.get_z()))


@property
def chain_str(self):
    return _char_array_to_str_list(self.get_chain())


@property
def name_str(self):
    return _char_array_to_str_list(self.get_name())


@property
def resname_str(self):
    return _char_array_to_str_list(self.get_resname())


@property
def elem_str(self):
    return _char_array_to_str_list(self.get_elem())


@property
def alterloc_str(self):
    return _char_array_to_str_list(self.get_alterloc())


@property
def insertres_str(self):
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
    return _FieldProxy(self, 'get_resname', 'set_resname')


@property
def resid(self):
    return _FieldProxy(self, 'get_resid', 'set_resid')


@property
def uniq_resid(self):
    return _FieldProxy(self, 'get_uniqresid', 'set_uniqresid')


@property
def chain(self):
    return _FieldProxy(self, 'get_chain', 'set_chain')


@property
def name(self):
    return _FieldProxy(self, 'get_name', 'set_name')


@property
def num(self):
    return _FieldProxy(self, 'get_num', 'set_num')


@property
def x(self):
    return _FieldProxy(self, 'get_x', 'set_x')


@property
def y(self):
    return _FieldProxy(self, 'get_y', 'set_y')


@property
def z(self):
    return _FieldProxy(self, 'get_z', 'set_z')


@property
def beta(self):
    return _FieldProxy(self, 'get_beta', 'set_beta')


@property
def occ(self):
    return _FieldProxy(self, 'get_occ', 'set_occ')


@property
def xyz(self):
    return np.column_stack((self.get_x(), self.get_y(), self.get_z()))


@property
def chain_str(self):
    return _char_array_to_str_list(self.get_chain())


@property
def name_str(self):
    return _char_array_to_str_list(self.get_name())


@property
def resname_str(self):
    return _char_array_to_str_list(self.get_resname())


@property
def elem_str(self):
    return _char_array_to_str_list(self.get_elem())


@property
def alterloc_str(self):
    return _char_array_to_str_list(self.get_alterloc())


@property
def insertres_str(self):
    return _char_array_to_str_list(self.get_insertres())


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
