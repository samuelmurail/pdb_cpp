[![Documentation Status](https://readthedocs.org/projects/pdb-cpp/badge/?version=latest)](https://pdb-cpp.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/samuelmurail/pdb_cpp/graph/badge.svg?token=QWFM35BX2N)](https://codecov.io/gh/samuelmurail/pdb_cpp)
[![CI](https://github.com/samuelmurail/pdb_cpp/actions/workflows/ci.yml/badge.svg)](https://github.com/samuelmurail/pdb_cpp/actions/workflows/ci.yml)
[![PyPI - Version](https://img.shields.io/pypi/v/pdb-cpp)](https://pypi.org/project/pdb-cpp/)
[![Downloads](https://static.pepy.tech/badge/pdb-cpp)](https://pepy.tech/project/pdb-cpp)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)


# pdb_cpp

`pdb_cpp` is a structural bioinformatics toolkit with a C++ core and Python API for fast PDB/mmCIF parsing, atom selection, sequence/structure alignment, TM-score, and DockQ evaluation.

<img src="https://raw.githubusercontent.com/samuelmurail/pdb_cpp/master/docs/source/logo.png" alt="pdb_cpp logo" width="300" style="display: block; margin: auto;"/>

## What is included

- Read/write `.pdb`, `.cif`, `.pqr`, and `.gro` files
- Atom/residue/chain selections (including geometric `within` queries)
- Sequence extraction and pairwise sequence alignment
- Sequence-based structural superposition and chain-permutation alignment
- TM-align/TM-score through the bundled USalign/TM-align core
- DockQ metrics (`DockQ`, `Fnat`, `Fnonnat`, `LRMS`, `iRMS`, `rRMS`)
- Hydrogen bond detection (Baker & Hubbard geometric method, no explicit H required)
- Solvent-accessible surface area (Shrake-Rupley) on `Model`, plus buried protein-protein surface and shape-complementarity helpers
- Secondary structure assignment
- Core geometric helpers (e.g., distance matrix)

## Installation

### From PyPI

```bash
python -m pip install pdb-cpp
```

### From source

```bash
git clone https://github.com/samuelmurail/pdb_cpp
cd pdb_cpp
python -m pip install -e .
```

For development:

```bash
python -m pip install -r requirements.txt
pytest
```

## Quick start

```python
from pdb_cpp import Coor

# Load from local file
coor = Coor("tests/input/1y0m.cif")

# Or fetch by PDB ID (mmCIF is downloaded and cached)
coor_pdb = Coor(pdb_id="1y0m")

# Or use the RCSB helper for explicit structure choices
from pdb_cpp import rcsb

bio_assembly = rcsb.load("5a9z", structure="biological_assembly", assembly_id=1)
asym_unit = rcsb.load("5a9z", structure="asymmetric_unit")

print(coor.model_num)        # number of models
print(coor.get_aa_seq())     # chain -> sequence

# Write selection/structure back to disk
coor.write("out_structure.pdb")
```

## Documentation map

For complete usage documentation, use the project docs site and source pages:

- Basic tutorial: `docs/source/basic_example.md`
- Full feature guide: `docs/source/functionality.md`
- Copy-paste recipes: `docs/source/quick_recipes.md`
- Installation and build notes: `docs/source/installation.md`
- API reference entry: `docs/source/pdb_cpp.rst`

Online docs: https://samuelmurail.github.io/pdb_cpp/

## Documentation

- API and examples: https://samuelmurail.github.io/pdb_cpp/
- Docs sources: `docs/source/`

## Notes for contributors (C++ core)

When adding C++ features:

1. Add implementation files in `src/pdb_cpp/_core/`
2. Register sources in `setup.py`
3. Expose bindings in `src/pdb_cpp/_core/pybind.cpp`
4. Reinstall extension (`pip install -e . --no-build-isolation`) and run tests