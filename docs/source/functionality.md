# Functionality Guide

This guide documents every major `pdb_cpp` workflow with fully worked
examples. For copy-paste snippets see [Quick Recipes](quick_recipes.md);
for the auto-generated API reference see [pdb_cpp package](pdb_cpp.rst).

---

## 1. Reading and writing structures

### Supported formats

| Format | Extension | Read | Write | Notes                                    |
|--------|-----------|------|-------|------------------------------------------|
| PDB    | `.pdb`    | yes  | yes   |                                          |
| mmCIF  | `.cif`    | yes  | yes   |                                          |
| PQR    | `.pqr`    | yes  | yes   | Charge in `occ`, radius in `beta` slot   |
| GRO    | `.gro`    | yes  | yes   | GROMACS format; coordinates converted nm ↔ Å |

### Loading from a local file

```python
from pdb_cpp import Coor

coor = Coor("structure.cif")   # auto-detects format by extension
coor = Coor("structure.pdb")
```

### Fetching from the PDB archive

```python
coor = Coor(pdb_id="1y0m")    # downloads mmCIF, caches in /tmp/pdb_cpp_cache/
```

The cached file is reused on subsequent calls with the same PDB ID.

### Writing output

```python
coor.write("output.pdb")      # PDB format
coor.write("output.cif")      # mmCIF format
coor.write("output.pqr")      # PQR format
coor.write("output.gro")      # GRO format
```

The format is determined by the file extension. Selections can also be
written directly:

```python
sel = coor.select_atoms("chain A")
sel.write("chain_A.pdb")
```

---

## 2. The `Coor` and `Model` objects

A `Coor` object holds one or more `Model` frames. Many properties are
available on both classes via Python property accessors.

### Key `Coor` properties

| Property          | Type            | Description                                   |
|-------------------|-----------------|-----------------------------------------------|
| `len`             | `int`           | Number of atoms in the active model           |
| `model_num`       | `int`           | Number of models                              |
| `models`          | `list[Model]`   | All models                                    |
| `active_model`    | `int`           | Index of the active model (read/write)        |
| `conect`          | `list`          | CONECT records (PDB format)                   |

### Per-atom properties (on both `Coor` and `Model`)

These properties return data for the active model when called on `Coor`,
or directly when called on a `Model`.

| Property        | Type            | Description                              |
|-----------------|-----------------|------------------------------------------|
| `x`, `y`, `z`  | `list[float]`   | Cartesian coordinates                    |
| `xyz`           | `ndarray(N,3)`  | Coordinates as a NumPy array             |
| `name`          | `list`          | Atom names (fixed-length char arrays)    |
| `name_str`      | `list[str]`     | Atom names as Python strings             |
| `resname`       | `list`          | Residue names (char arrays)              |
| `resname_str`   | `list[str]`     | Residue names as strings                 |
| `resid`         | `list[int]`     | Residue sequence numbers                 |
| `uniq_resid`    | `list[int]`     | Unique (0-based) residue identifiers     |
| `chain`         | `list`          | Chain identifiers (char arrays)          |
| `chain_str`     | `list[str]`     | Chain identifiers as strings             |
| `num`           | `list[int]`     | Atom serial numbers                      |
| `beta`          | `list[float]`   | B-factors (temperature factors)          |
| `occ`           | `list[float]`   | Occupancies                              |
| `elem_str`      | `list[str]`     | Element symbols                          |
| `alterloc_str`  | `list[str]`     | Alternate location indicators            |
| `insertres_str` | `list[str]`     | Residue insertion codes                  |

```python
coor = Coor("tests/input/1y0m.cif")

# Access via Coor (delegates to active model)
print(coor.x[:5])
print(coor.chain_str[:5])
print(coor.xyz.shape)

# Access a specific model directly
model_0 = coor.models[0]
print(model_0.name_str[:5])
print(model_0.beta[:5])
```

### Multi-model structures

```python
coor = Coor("tests/input/2rri.pdb")
print(coor.model_num)   # e.g. 20 models (NMR ensemble)

# Switch active model
coor.active_model = 5
print(coor.x[:3])       # coordinates from model 5

# Iterate over all models
for i, model in enumerate(coor.models):
    print(f"Model {i}: {model.len} atoms")
```

### Centroid calculation

```python
model = coor.models[0]

# Centroid of all atoms
centroid = model.get_centroid()

# Centroid of specific atom indices
ca_indices = coor.get_index_select("name CA")
centroid_ca = model.get_centroid(ca_indices)
```

---

## 3. Atom selection

### Selection syntax

Selections are performed with `coor.select_atoms(selection_string)` which
returns a new `Coor` object containing only the matching atoms.

#### Field selectors

| Keyword    | Matches                          | Example                          |
|------------|----------------------------------|----------------------------------|
| `name`     | Atom name                        | `name CA`                        |
| `resname`  | Residue name                     | `resname ALA GLY`                |
| `chain`    | Chain identifier                 | `chain A B`                      |
| `resid`    | Residue sequence number          | `resid 1 2 3`                    |
| `residue`  | Unique residue index (0-based)   | `residue >= 10 and residue < 50` |
| `altloc`   | Alternate location indicator     | `not altloc B C D`               |
| `x`, `y`, `z` | Coordinate value             | `x >= 20.0`                      |
| `beta`     | B-factor                         | `beta < 50.0`                    |
| `occ`      | Occupancy                        | `occ > 0.5`                      |

#### Shortcut keywords

| Keyword    | Equivalent to                               |
|------------|---------------------------------------------|
| `protein`  | Standard amino acid residue names            |
| `backbone` | `name CA C N O` within protein residues      |
| `noh`      | Non-hydrogen atoms                           |

#### Logical operators

- `and` — intersection
- `or` — union
- `not` — negation

#### Comparison operators

`==`, `!=`, `<`, `<=`, `>`, `>=` — used with numeric fields (`resid`,
`residue`, `x`, `y`, `z`, `beta`, `occ`).

#### Spatial queries

```
within <distance> of <sub-selection>
```

Selects atoms within `distance` Ångströms of any atom matching the
sub-selection.

### Selection examples

```python
from pdb_cpp import Coor

coor = Coor("tests/input/1rxz.pdb")

# Chain selection
chain_a = coor.select_atoms("chain A")

# Backbone of a residue range
seg = coor.select_atoms("backbone and chain A and residue >= 6 and residue <= 58")

# Interface: atoms of A within 5 Å of B
interface = coor.select_atoms("chain A and within 5.0 of chain B")

# Exclude water and alternate locations
clean = coor.select_atoms("protein and not altloc B C D and not resname HOH")

# Numeric coordinate filter
slab = coor.select_atoms("name CA and x >= 10.0 and x <= 30.0")
```

### Getting indices without creating a new `Coor`

```python
indices = coor.get_index_select("name CA and chain A")
# Returns list[int] of atom indices in the active model
```

### Frame-specific selections (multi-model)

```python
# Select using coordinates from a specific model frame
sel = coor.select_atoms("within 5.0 of chain A", frame=3)
```

---

## 4. Sequence extraction

### Amino acid sequences

```python
seqs = coor.get_aa_seq()
# Returns dict[str, str]: chain ID -> one-letter sequence
print(seqs["A"])
```

With gap insertion for missing residues:

```python
seqs = coor.get_aa_seq(gap_in_seq=True)   # default: gaps inserted
seqs = coor.get_aa_seq(gap_in_seq=False)   # no gap characters
```

### D/L amino acid sequences

```python
dl_seqs = coor.get_aa_DL_seq()
# D-amino acids are returned as lowercase letters
```

### Amino acid and nucleic acid sequences

```python
all_seqs = coor.get_aa_na_seq()
# Includes both protein and nucleic acid chains
```

### Via the `sequence` module

```python
from pdb_cpp import sequence

seqs = sequence.get_aa_seq(coor)
dl_seqs = sequence.get_aa_DL_seq(coor)
```

---

## 5. Sequence alignment

`pdb_cpp.alignment` provides pairwise sequence alignment using a BLOSUM62
substitution matrix and affine gap penalties.

```python
from pdb_cpp import alignment

seq_1 = "ACDEFGHIKLMNPQRSTVWY"
seq_2 = "ACDFGHIKLMNPQRSTVWY"   # missing E

aln_1, aln_2, score = alignment.align_seq(seq_1, seq_2)
print(aln_1)
print(aln_2)
print(f"Score: {score}")

# Pretty-print with identity/similarity stats
alignment.print_align_seq(aln_1, aln_2)
```

### Parameters

| Parameter     | Default | Description                        |
|---------------|---------|------------------------------------|
| `gap_cost`    | -11     | Gap opening penalty                |
| `gap_ext`     | -1      | Gap extension penalty              |
| `matrix_file` | `None`  | Custom substitution matrix file; defaults to bundled BLOSUM62 |

---

## 6. Structural alignment (RMSD-based)

### Step-by-step alignment

```python
from pdb_cpp import Coor, core, analysis

coor_1 = Coor("tests/input/1u85.pdb")
coor_2 = Coor("tests/input/1ubd.pdb")

# 1. Find common atoms by sequence alignment
idx_1, idx_2 = core.get_common_atoms(
    coor_1, coor_2,
    chain_1=["A"], chain_2=["C"],
    back_names=["C", "N", "O", "CA"],  # backbone atoms
)

# 2. Superpose coor_1 onto coor_2 (modifies coor_1 in-place)
core.coor_align(coor_1, coor_2, idx_1, idx_2, frame_ref=0)

# 3. Compute RMSD
rmsd_values = analysis.rmsd(coor_1, coor_2, index_list=[idx_1, idx_2])
print(f"RMSD: {rmsd_values[0]:.3f} Å")
```

### One-step sequence-based alignment

```python
rmsd_list, align_idx_1, align_idx_2 = core.align_seq_based(
    coor_1, coor_2,
    chain_1=["A"], chain_2=["C"],
)
print(f"RMSD: {rmsd_list[0]:.3f} Å")
```

This combines `get_common_atoms`, `coor_align`, and RMSD in a single call.

### Chain-permutation alignment (multi-chain complexes)

When chain correspondence is unknown:

```python
from pdb_cpp import alignment

rmsds, mappings = alignment.align_chain_permutation(coor_1, coor_2)
print(f"Best RMSD: {rmsds[0]:.3f} Å")
```

This tries all possible chain pairings and returns the best RMSD.

### RMSD functions

```python
from pdb_cpp import analysis

# RMSD using a selection string
rmsd_values = analysis.rmsd(coor_1, coor_2, selection="name CA")

# RMSD with explicit index lists (no re-selection)
rmsd_values = analysis.rmsd(coor_1, coor_2, index_list=[idx_1, idx_2])
```

---

## 7. TM-align and TM-score

TM-align provides length-independent structural similarity scoring via
the bundled USalign implementation.

```python
from pdb_cpp import Coor
from pdb_cpp.core import tmalign_ca

coor_1 = Coor("tests/input/1y0m.cif")
coor_2 = Coor("tests/input/1ubd.pdb")

result = tmalign_ca(
    coor_1, coor_2,
    chain_1=["A"], chain_2=["C"],
    mm=1,    # 1 = TM-align mode
)
```

### `TMalignResult` fields

| Field    | Type    | Description                                         |
|----------|---------|-----------------------------------------------------|
| `TM1`    | `float` | TM-score normalized by length of structure 1        |
| `TM2`    | `float` | TM-score normalized by length of structure 2        |
| `TM_ali` | `float` | TM-score normalized by aligned length               |
| `rmsd`   | `float` | RMSD of aligned residues                             |
| `L_ali`  | `int`   | Number of aligned residues                           |
| `Liden`  | `int`   | Number of identical aligned residues                 |
| `seqxA`  | `str`   | Aligned sequence of structure 1                      |
| `seqyA`  | `str`   | Aligned sequence of structure 2                      |
| `seqM`   | `str`   | Alignment match string (`:` aligned, `.` close, ` ` far) |

### The `mm` parameter

| Value | Meaning                  |
|-------|--------------------------|
| 0     | Monomer TM-align         |
| 1     | Monomer TM-align (same)  |

### Multi-chain alignment

For multi-chain alignment with chain-pairing, pass multiple chains:

```python
result = tmalign_ca(
    coor_1, coor_2,
    chain_1=["A", "B"], chain_2=["A", "B"],
    mm=1,
)
```

### Citation

If you use the TM-align functionality, please cite:

- Chengxin Zhang, Morgan Shine, Anna Marie Pyle, Yang Zhang (2022) *Nat Methods*. 19(9), 1109-1115.
- Chengxin Zhang, Anna Marie Pyle (2022) *iScience*. 25(10), 105218.

---

## 8. Secondary structure assignment

Uses the TM-align secondary structure assignment algorithm (DSSP-like).

```python
from pdb_cpp import Coor, TMalign

coor = Coor("tests/input/1y0m.cif")
ss_list = TMalign.compute_secondary_structure(coor)
```

Returns a list (one entry per model) of dictionaries mapping chain ID to
a secondary structure string.

### Secondary structure codes

| Code | Meaning          |
|------|------------------|
| `H`  | Helix            |
| `E`  | Strand (sheet)   |
| `T`  | Turn             |
| `C`  | Coil / other     |

```python
for chain_id, ss_string in ss_list[0].items():
    print(f"Chain {chain_id}: {ss_string}")
```

---

## 9. DockQ scoring

DockQ evaluates docking model quality by combining interface contacts,
ligand RMSD, and interface RMSD into a single score.

### Basic usage (automatic chain-role inference)

```python
from pdb_cpp import Coor, analysis

model = Coor("tests/input/1rxz_colabfold_model_1.pdb")
native = Coor("tests/input/1rxz.pdb")

scores = analysis.dockQ(model, native)
```

### With explicit chain roles

```python
scores = analysis.dockQ(
    model, native,
    rec_chains=["B"],
    lig_chains=["C"],
    native_rec_chains=["A"],
    native_lig_chains=["B"],
)
```

### DockQ output dictionary

| Key        | Type          | Description                          |
|------------|---------------|--------------------------------------|
| `DockQ`    | `list[float]` | Combined DockQ score (0-1)           |
| `Fnat`     | `list[float]` | Fraction of native contacts          |
| `Fnonnat`  | `list[float]` | Fraction of non-native contacts      |
| `LRMS`     | `list[float]` | Ligand RMSD (Å)                      |
| `iRMS`     | `list[float]` | Interface RMSD (Å)                   |
| `rRMS`     | `list[float]` | Receptor RMSD (Å)                    |

Each value is a list with one entry per model in the input `Coor`.

### DockQ quality thresholds

| DockQ range | Quality    |
|-------------|------------|
| 0.00 – 0.23 | Incorrect  |
| 0.23 – 0.49 | Acceptable |
| 0.49 – 0.80 | Medium     |
| 0.80 – 1.00 | High       |

### Individual metric functions

```python
# Fraction of native/non-native contacts
fnat_list, fnonnat_list = analysis.native_contact(
    model, native,
    rec_chains=["A"], lig_chains=["B"],
    native_rec_chains=["A"], native_lig_chains=["B"],
    cutoff=5.0,
)

# Interface RMSD
irmsd_list = analysis.interface_rmsd(
    model, native,
    rec_chains_native=["A"],
    lig_chains_native=["B"],
    cutoff=10.0,
)
```

### Citation

If you use DockQ scoring, please cite:

- Mirabello, C. & Wallner, B. (2024) *Bioinformatics*. DOI: 10.1093/bioinformatics/btae586

---

## 10. Geometry utilities

### Distance matrix

```python
from pdb_cpp import Coor, geom

coor = Coor("tests/input/1y0m.cif")
ca = coor.select_atoms("name CA")

# Pairwise distance matrix between two coordinate sets
dmat = geom.distance_matrix(ca, ca)
print(dmat.shape)  # (N_ca, N_ca)
```

The function accepts any two `Coor` objects (or selections) and returns an
`(N, M)` float32 NumPy array of Euclidean distances.

### Hybrid-36 encoding/decoding

PDB format uses hybrid-36 encoding for atom/residue numbers > 99999:

```python
from pdb_cpp.core import hy36encode, hy36decode

encoded = hy36encode(5, 100000)   # "A0000"
decoded = hy36decode(5, "A0000")  # 100000
```

---

## 11. Data cleaning utilities

### Remove incomplete backbone residues

Residues missing any backbone atom (CA, C, N, O by default) are removed:

```python
from pdb_cpp import select

clean_coor = select.remove_incomplete_backbone_residues(coor)

# Or with custom backbone definition
clean_coor = select.remove_incomplete_backbone_residues(
    coor, back_atom=["CA", "C", "N"]
)
```

This is recommended before structural alignment to avoid mismatches.

---

## 12. Module summary

| Module              | Key functions / classes                                   |
|---------------------|-----------------------------------------------------------|
| `pdb_cpp`           | `Coor`, `Model` — main data objects                           |
| `pdb_cpp.core`      | `get_common_atoms`, `coor_align`, `align_seq_based`, `tmalign_ca`, `distance_matrix`, `compute_SS`, `hy36encode`, `hy36decode` |
| `pdb_cpp.alignment` | `align_seq`, `print_align_seq`, `align_chain_permutation` |
| `pdb_cpp.analysis`  | `rmsd`, `interface_rmsd`, `native_contact`, `dockQ`       |
| `pdb_cpp.TMalign`   | `compute_secondary_structure`                             |
| `pdb_cpp.geom`      | `distance_matrix`                                         |
| `pdb_cpp.select`    | `remove_incomplete_backbone_residues`                     |
| `pdb_cpp.sequence`  | `get_aa_seq`, `get_aa_DL_seq`                             |
