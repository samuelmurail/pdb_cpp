#!/usr/bin/env python3
# coding: utf-8

"""High-level analysis namespace.

This package groups structure-analysis helpers into topic-oriented modules:

- :mod:`pdb_cpp.analysis.dockq`
- :mod:`pdb_cpp.analysis.sasa`
- :mod:`pdb_cpp.analysis.hbonds`

The historical flat API is preserved, so existing code such as
``pdb_cpp.analysis.rmsd(...)`` keeps working.
"""

from . import dockq, interaction, salt_bridge as _salt_bridge_module
from . import hbonds as _hbonds_module
from . import sasa as _sasa_module
from .dockq import dockQ, dockQ_multimer, interface_rmsd, native_contact, rmsd
from .hbonds import hbonds as compute_hbonds
from .salt_bridge import salt_bridges
from .sasa import buried_surface_area


class _CallableNamespace:
    def __init__(self, func, module):
        self._func = func
        self._module = module

    def __call__(self, *args, **kwargs):
        return self._func(*args, **kwargs)

    def __getattr__(self, name):
        return getattr(self._module, name)


sasa = _CallableNamespace(_sasa_module.sasa, _sasa_module)
hbonds = _CallableNamespace(_hbonds_module.hbonds, _hbonds_module)
salt_bridge = _salt_bridge_module

__all__ = [
    "dockq",
    "interaction",
    "sasa",
    "hbonds",
    "salt_bridge",
    "rmsd",
    "interface_rmsd",
    "native_contact",
    "dockQ",
    "dockQ_multimer",
    "buried_surface_area",
    "compute_hbonds",
    "salt_bridges",
]