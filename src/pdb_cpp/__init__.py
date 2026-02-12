#!/usr/bin/env python3
# coding: utf-8

# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2025, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__version__ = "0.0.1"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Beta"

import os
import tempfile
import urllib.request
from urllib.error import HTTPError

import numpy as np

from .core import Coor, Model

_coor_init = Coor.__init__


def _fetch_mmcif(pdb_id):
    pdb_id = pdb_id.lower()
    cache_dir = os.path.join(tempfile.gettempdir(), "pdb_cpp_cache")
    os.makedirs(cache_dir, exist_ok=True)
    local_path = os.path.join(cache_dir, f"{pdb_id}.cif")
    if os.path.exists(local_path):
        return local_path

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    try:
        with urllib.request.urlopen(url) as response:
            data = response.read()
    except HTTPError as exc:
        raise ValueError(f"Failed to fetch mmCIF for PDB ID {pdb_id}") from exc

    with open(local_path, "wb") as handle:
        handle.write(data)
    return local_path


def _coor_init_with_pdb_id(self, coor_in=None, pdb_id=None):
    _coor_init(self)
    if coor_in is not None and pdb_id is not None:
        raise ValueError("Provide only one of coor_in or pdb_id")
    if pdb_id is not None:
        local_path = _fetch_mmcif(pdb_id)
        self.read(local_path)
    elif coor_in is not None:
        self.read(coor_in)


Coor.__init__ = _coor_init_with_pdb_id


def _char_array_to_str_list(array_like):
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
    """Return the atom number of the selection."""
    return self.get_num()


@property
def residue(self):
    """Return the residue of the selection."""
    return self.get_uniqresid()


@property
def uniq_resid(self):
    """Return the residue of the selection."""
    return self.get_uniqresid()


@property
def name(self):
    """Return the atom name of the selection."""
    return self.get_name()


@property
def chain(self):
    """Return the chain of the selection."""
    return self.get_chain()


@property
def resname(self):
    """Return the residue name of the selection."""
    return self.get_resname()


@property
def resid(self):
    """Return the residue id of the selection."""
    return self.get_resid()


@property
def x(self):
    """Return the x coordinate of the selection."""
    return self.get_x()


@property
def y(self):
    """Return the y coordinate of the selection."""
    return self.get_y()


@property
def z(self):
    """Return the z coordinate of the selection."""
    return self.get_z()


@property
def beta(self):
    """Return the beta factor of the selection."""
    return self.get_beta()


@property
def occ(self):
    """Return the occupancy of the selection."""
    return self.get_occupancy()


@property
def xyz(self):
    """Return the xyz coordinates as an (N, 3) array."""
    return np.column_stack((self.get_x(), self.get_y(), self.get_z()))


@property
def chain_str(self):
    """Return the chain identifiers as strings."""
    return _char_array_to_str_list(self.get_chain())


@property
def name_str(self):
    """Return the atom names as strings."""
    return _char_array_to_str_list(self.get_name())


@property
def resname_str(self):
    """Return the residue names as strings."""
    return _char_array_to_str_list(self.get_resname())


@property
def elem_str(self):
    """Return the element symbols as strings."""
    return _char_array_to_str_list(self.get_elem())


@property
def alterloc_str(self):
    """Return the alternate location identifiers as strings."""
    return _char_array_to_str_list(self.get_alterloc())


@property
def insertres_str(self):
    """Return the residue insertion codes as strings."""
    return _char_array_to_str_list(self.get_insertres())


Model.num = num
Model.residue = residue
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
    """Return the number of atoms in the selection."""
    return self.size()


Model.len = length
Coor.len = length


@property
def models(self):
    """Return the model of the selection."""
    return self.get_all_Models()


@property
def model_num(self):
    """Return the model of the selection."""
    return self.model_size()


# @property
# def residue(self):
#     """Return the residue of the selection."""
#     return self.models[0].get_uniqresid()

Coor.models = models
Coor.model_num = model_num


@property
def resname(self):
    """Return the residue name of the selection."""
    return self.models[self.active_model].get_resname()


@property
def resid(self):
    """Return the residue id of the selection."""
    return self.models[self.active_model].get_resid()


@property
def residue(self):
    """Return the residue of the selection."""
    return self.models[self.active_model].get_uniqresid()


@property
def uniq_resid(self):
    """Return the residue of the selection."""
    return self.models[self.active_model].get_uniqresid()


@property
def chain(self):
    """Return the chain of the selection."""
    return self.models[self.active_model].get_chain()


@property
def name(self):
    """Return the atom name of the selection."""
    return self.models[self.active_model].get_name()


@property
def num(self):
    """Return the atom number of the selection."""
    return self.models[self.active_model].get_num()


@property
def x(self):
    """Return the x coordinate of the selection."""
    return self.models[self.active_model].get_x()


@property
def y(self):
    """Return the y coordinate of the selection."""
    return self.models[self.active_model].get_y()


@property
def z(self):
    """Return the z coordinate of the selection."""
    return self.models[self.active_model].get_z()


@property
def beta(self):
    """Return the beta factor of the selection."""
    return self.models[self.active_model].get_beta()


@property
def occ(self):
    """Return the occupancy of the selection."""
    return self.models[self.active_model].get_occupancy()


@property
def xyz(self):
    """Return the xyz coordinates as an (N, 3) array."""
    model = self.models[self.active_model]
    return np.column_stack((model.get_x(), model.get_y(), model.get_z()))


@property
def chain_str(self):
    """Return the chain identifiers as strings."""
    return _char_array_to_str_list(self.models[self.active_model].get_chain())


@property
def name_str(self):
    """Return the atom names as strings."""
    return _char_array_to_str_list(self.models[self.active_model].get_name())


@property
def resname_str(self):
    """Return the residue names as strings."""
    return _char_array_to_str_list(self.models[self.active_model].get_resname())


@property
def elem_str(self):
    """Return the element symbols as strings."""
    return _char_array_to_str_list(self.models[self.active_model].get_elem())


@property
def alterloc_str(self):
    """Return the alternate location identifiers as strings."""
    return _char_array_to_str_list(self.models[self.active_model].get_alterloc())


@property
def insertres_str(self):
    """Return the residue insertion codes as strings."""
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
Coor.residue = residue
Coor.uniq_resid = uniq_resid

from . import sequence
