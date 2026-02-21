pdb\_cpp package
================

The package is organized around a high-level ``Coor`` object and focused helper
modules.

Main workflow modules:

- ``pdb_cpp.core``: C++-accelerated core operations (I/O, selection, atom mapping,
  coordinate alignment, TM-align, sequence-based structural alignment).
- ``pdb_cpp.analysis``: RMSD, interface RMSD, native contacts, DockQ.
- ``pdb_cpp.alignment``: sequence alignment and chain-permutation helpers.
- ``pdb_cpp.TMalign``: secondary structure helper built on TM-align core.
- ``pdb_cpp.sequence``: sequence extraction wrappers.
- ``pdb_cpp.geom``: geometry helpers (distance matrix).

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   pdb_cpp.data

Submodules
----------

pdb\_cpp.TMalign module
-----------------------

.. automodule:: pdb_cpp.TMalign
   :members:
   :show-inheritance:
   :undoc-members:

pdb\_cpp.alignment module
-------------------------

.. automodule:: pdb_cpp.alignment
   :members:
   :show-inheritance:
   :undoc-members:

pdb\_cpp.analysis module
------------------------

.. automodule:: pdb_cpp.analysis
   :members:
   :show-inheritance:
   :undoc-members:

pdb\_cpp.core extension module
------------------------------

.. automodule:: pdb_cpp.core
   :members:
   :show-inheritance:
   :undoc-members:

pdb\_cpp.geom module
--------------------

.. automodule:: pdb_cpp.geom
   :members:
   :show-inheritance:
   :undoc-members:

pdb\_cpp.select module
----------------------

.. automodule:: pdb_cpp.select
   :members:
   :show-inheritance:
   :undoc-members:

pdb\_cpp.sequence module
------------------------

.. automodule:: pdb_cpp.sequence
   :members:
   :show-inheritance:
   :undoc-members:

Module contents
---------------

.. automodule:: pdb_cpp
   :members:
   :show-inheritance:
   :undoc-members:
