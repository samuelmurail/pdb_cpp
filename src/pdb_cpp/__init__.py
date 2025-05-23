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

from .core import Coor, Model

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

@property
def len(self):
    """Return the number of atoms in the selection."""
    return self.size()

Model.len = len
Coor.len = len

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
Coor.residue = residue
Coor.uniq_resid = uniq_resid