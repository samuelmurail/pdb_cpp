# Basic Usage

## Import and load structures

```python
from pdb_cpp import Coor

# Load from local file (PDB, mmCIF, PQR, or GRO — auto-detected by extension)
coor = Coor("tests/input/1y0m.cif")

# Load from PDB ID (downloaded as mmCIF and cached locally)
coor_from_id = Coor(pdb_id="1y0m")
```

`Coor` can read `.pdb`, `.cif`, `.pqr`, and `.gro` files and can contain one or more models.

## Read basic properties

```python
print(coor.model_num)            # number of models
print(coor.len)                  # number of atoms in active model
print(coor.get_uniq_chain())     # list of chain IDs
print(coor.get_aa_seq())         # dict: chain -> sequence
```

## Access atom-level data

Properties are available directly on both `Coor` (delegates to the active
model) and on individual `Model` objects:

```python
# Coordinates as NumPy array
print(coor.xyz[:5])              # shape (N, 3)

# Atom names, residue names, chains as Python strings
print(coor.name_str[:5])         # ['N', 'CA', 'C', 'O', 'CB']
print(coor.resname_str[:5])      # ['THR', 'THR', 'THR', 'THR', 'THR']
print(coor.chain_str[:5])        # ['A', 'A', 'A', 'A', 'A']

# Numeric properties
print(coor.resid[:5])            # residue sequence numbers
print(coor.beta[:5])             # B-factors
print(coor.occ[:5])              # occupancies

# Access a specific model
model = coor.models[0]
print(model.x[:3])
print(model.elem_str[:3])
```

## Working with multi-model structures

```python
coor = Coor("tests/input/2rri.pdb")
print(f"{coor.model_num} models, {coor.len} atoms each")

# Switch the active model
coor.active_model = 5
print(coor.xyz[:3])              # coordinates from model 5

# Iterate over all models
for i, model in enumerate(coor.models):
    print(f"Model {i}: {model.len} atoms, centroid = {model.get_centroid()}")
```

## Atom and residue selections

```python
# Chain selection
chain_a = coor.select_atoms("chain A")

# Standard protein/backbone keywords
protein_backbone = coor.select_atoms("protein and backbone")

# Residue ranges
segment = coor.select_atoms("chain A and residue >= 6 and residue <= 58")

# Spatial query: atoms in chain A within 5 Å of chain B
interface = coor.select_atoms("chain A and within 5.0 of chain B")

# Combine logic and numeric filters
filtered = coor.select_atoms("name CA and not within 5.0 of resname HOH and x >= 20.0")

# Get atom indices without creating a new Coor
indices = coor.get_index_select("name CA and chain A")
```

Supported selection primitives include:

- atom fields: `name`, `altloc`, `resname`, `chain`, `resid`, `residue`, `x`, `y`, `z`, `occ`, `beta`
- shortcuts: `protein`, `backbone`, `noh`
- operators: `and`, `or`, `not`, `within`, `==`, `!=`, `>`, `>=`, `<`, `<=`

## Write output structures

```python
interface.write("interface_chainA.pdb")
segment.write("segment_chainA.cif")
```

## Sequence alignment and coordinate superposition

```python
from pdb_cpp import alignment, core, analysis

coor_1 = Coor("tests/input/1u85.pdb")
coor_2 = Coor("tests/input/1ubd.pdb")

# Pairwise sequence alignment
seq_1 = coor_1.get_aa_seq()["A"]
seq_2 = coor_2.get_aa_seq()["C"]
aln_1, aln_2, score = alignment.align_seq(seq_1, seq_2)
alignment.print_align_seq(aln_1, aln_2)

# Find common backbone atoms and superpose
idx_1, idx_2 = core.get_common_atoms(coor_1, coor_2, chain_1=["A"], chain_2=["C"])
core.coor_align(coor_1, coor_2, idx_1, idx_2, frame_ref=0)

# Compute RMSD after alignment
rmsd = analysis.rmsd(coor_1, coor_2, index_list=[idx_1, idx_2])
print(f"RMSD: {rmsd[0]:.3f} Å")

# Or do everything in one step
rmsd_list, _, _ = core.align_seq_based(coor_1, coor_2, chain_1=["A"], chain_2=["C"])
```

## TM-score / TM-align

```python
from pdb_cpp.core import tmalign_ca

coor_1 = Coor("tests/input/1y0m.cif")
coor_2 = Coor("tests/input/1ubd.pdb")

tm = tmalign_ca(coor_1, coor_2, chain_1=["A"], chain_2=["C"], mm=1)

print(f"Aligned: {tm.L_ali} residues, RMSD: {tm.rmsd:.3f}")
print(f"TM-score (norm by 1): {tm.TM1:.4f}")
print(f"TM-score (norm by 2): {tm.TM2:.4f}")
```

## DockQ scoring

```python
from pdb_cpp import analysis

model = Coor("tests/input/1rxz_colabfold_model_1.pdb")
native = Coor("tests/input/1rxz.pdb")

scores = analysis.dockQ(model, native)
print(f"DockQ: {scores['DockQ'][0]:.3f}")
print(f"Fnat: {scores['Fnat'][0]:.3f}, Fnonnat: {scores['Fnonnat'][0]:.3f}")
print(f"LRMS: {scores['LRMS'][0]:.3f}, iRMS: {scores['iRMS'][0]:.3f}")
```

## Secondary structure

```python
from pdb_cpp import TMalign

ss = TMalign.compute_secondary_structure(coor)
for chain_id, ss_string in ss[0].items():
    print(f"Chain {chain_id}: {ss_string}")
```

## Distance matrix

```python
from pdb_cpp import geom

ca = coor.select_atoms("name CA")
dmat = geom.distance_matrix(ca, ca)
print(f"Distance matrix shape: {dmat.shape}")
```

