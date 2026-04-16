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

## Selection language (including complex selections)

`select_atoms()` supports boolean logic, numeric comparisons, and spatial queries:

```python
# Backbone residues 6..58 from chain A
sel_1 = coor.select_atoms("backbone and chain A and residue >= 6 and residue <= 58")

# Interface-like query: atoms of chain A within 5 Å of chain B
sel_2 = coor.select_atoms("chain A and within 5.0 of chain B")

# Combination with negation
sel_3 = coor.select_atoms("name CA and not within 5.0 of resname HOH")

# Numeric filters
sel_4 = coor.select_atoms("resname HOH and x >= 20.0")
```

Common keywords: `name`, `resname`, `chain`, `resid`, `residue`, `x`, `y`, `z`, `beta`, `occ`, `protein`, `backbone`, `noh`, `within`.

## Sequence alignment and structure superposition

```python
from pdb_cpp import Coor, alignment, core, analysis

coor_1 = Coor("tests/input/1u85.pdb")
coor_2 = Coor("tests/input/1ubd.pdb")

seq_1 = coor_1.get_aa_seq()["A"]
seq_2 = coor_2.get_aa_seq()["C"]
aln_1, aln_2, score = alignment.align_seq(seq_1, seq_2)
alignment.print_align_seq(aln_1, aln_2)

# Get atom correspondences and align coordinates in-place
idx_1, idx_2 = core.get_common_atoms(coor_1, coor_2, chain_1=["A"], chain_2=["C"])
core.coor_align(coor_1, coor_2, idx_1, idx_2, frame_ref=0)

# RMSD after alignment
rmsd_values = analysis.rmsd(coor_1, coor_2, index_list=[idx_1, idx_2])
print(rmsd_values[0])
```

For multi-chain complexes with uncertain chain mapping, use chain permutation:

```python
rmsds, index_mappings = alignment.align_chain_permutation(coor_1, coor_2)
```

## TM-align / TM-score

```python
from pdb_cpp import Coor
from pdb_cpp.core import tmalign_ca

ref = Coor("tests/input/1y0m.cif")
mob = Coor("tests/input/1ubd.pdb")

result = tmalign_ca(ref, mob, chain_1=["A"], chain_2=["C"], mm=1)

print(result.L_ali)  # aligned length
print(result.rmsd)   # RMSD on aligned residues
print(result.TM1)    # TM-score normalized by structure 1
print(result.TM2)    # TM-score normalized by structure 2
```

If you use the USalign/TM-align functionality in `pdb_cpp`, please cite:

- Chengxin Zhang, Morgan Shine, Anna Marie Pyle, Yang Zhang (2022) Nat Methods. 19(9), 1109-1115.
- Chengxin Zhang, Anna Marie Pyle (2022) iScience. 25(10), 105218.

Secondary structure helper:

```python
from pdb_cpp import TMalign

ss = TMalign.compute_secondary_structure(ref)
print(ss[0]["A"])  # DSSP-like secondary structure string for chain A
```

## DockQ scoring

```python
from pdb_cpp import Coor, analysis

model = Coor("tests/input/1rxz_colabfold_model_1.pdb")
native = Coor("tests/input/1rxz.pdb")

# Chain roles can be inferred automatically,
# or provided explicitly with rec_chains/lig_chains arguments.
scores = analysis.dockQ(model, native)

print(scores["DockQ"][0])
print(scores["Fnat"][0], scores["Fnonnat"][0])
print(scores["LRMS"][0], scores["iRMS"][0], scores["rRMS"][0])
```

Command-line usage is available through the installed `pdb_cpp_dockq` bin,
which uses `analysis.dockQ_multimer()` and therefore supports the 
automatic native-to-model chain mapping:

```bash
pdb_cpp_dockq tests/input/1a2k_model.pdb tests/input/1a2k.pdb

# Optional explicit native:model mapping
pdb_cpp_dockq tests/input/1a2k_model.pdb tests/input/1a2k.pdb --chain-map A:B,B:A,C:C
```

If you use DockQ scoring in `pdb_cpp`, please cite:

- DockQ, DOI: 10.1093/bioinformatics/btae586

## Hydrogen bond detection

`pdb_cpp.hbond` identifies hydrogen bonds using the Baker & Hubbard
geometric criteria. Hydrogen positions are reconstructed algebraically when
not present in the file, so no pre-processing step is required.

```python
from pdb_cpp import Coor
from pdb_cpp import hbond

coor = Coor("tests/input/2rri.cif")

# One list of HBond objects per model frame
all_bonds = hbond.hbonds(coor)
print(f"Model 0: {len(all_bonds[0])} H-bonds")

# Inspect bond geometry
b = all_bonds[0][0]
print(f"Donor  : {b.donor_chain}{b.donor_resid} {b.donor_heavy_name}")
print(f"Acceptor: {b.acceptor_chain}{b.acceptor_resid} {b.acceptor_name}")
print(f"d(D..A) = {b.dist_DA:.2f} Å  ∠DHA = {b.angle_DHA:.1f}°")

# Cross-selection: protein donors to nucleic-acid acceptors
rna_bonds = hbond.hbonds(coor, donor_sel="protein", acceptor_sel="nucleic")
```

Default cutoffs follow Baker & Hubbard (1984):
`dist_DA_cutoff=3.5 Å`, `dist_HA_cutoff=2.5 Å`, `angle_cutoff=90°`.

## Geometry utilities

```python
from pdb_cpp import Coor, geom

coor = Coor("tests/input/1y0m.cif")
ca = coor.select_atoms("name CA")
dist = geom.distance_matrix(ca, ca)
print(dist.shape)
```

## Documentation

- API and examples: https://samuelmurail.github.io/pdb_cpp/
- Docs sources: `docs/source/`

## Notes for contributors (C++ core)

When adding C++ features:

1. Add implementation files in `src/pdb_cpp/_core/`
2. Register sources in `setup.py`
3. Expose bindings in `src/pdb_cpp/_core/pybind.cpp`
4. Reinstall extension (`pip install -e . --no-build-isolation`) and run tests