# Functionality Guide

This guide documents every major `pdb_cpp` workflow with fully worked
examples. For copy-paste snippets see [Quick Recipes](quick_recipes.md);
for the auto-generated API reference see [pdb_cpp package](pdb_cpp.rst).

---

## 1. Reading and writing structures

### Supported formats

| Format | Extension | Read | Write | Notes                                    |
|--------|-----------|------|-------|------------------------------------------|
| PDB    | `.pdb`    | yes  | yes   | Bond topology via `CONECT` lines         |
| mmCIF  | `.cif`    | yes  | yes   | Bond topology via `_struct_conn` loop    |
| PQR    | `.pqr`    | yes  | yes   | Charge in `occ`, radius in `beta` slot   |
| GRO    | `.gro`    | yes  | yes   | GROMACS format; coordinates converted nm ↔ Å |

### Loading from a local file

```python
from pdb_cpp import Coor

coor = Coor("structure.cif")   # auto-detects format by extension
coor = Coor("structure.pdb")
```

The format can also be forced explicitly with the `format` keyword.
This is useful when the file extension is absent or misleading — for
example a compressed pipeline that writes to a generic filename:

```python
from pdb_cpp import Coor

# Via the constructor
coor = Coor("structure.dat", format="pdb")   # treat as PDB
coor = Coor("structure.dat", format="cif")   # treat as mmCIF
coor = Coor("structure.dat", format="pqr")   # treat as PQR
coor = Coor("structure.dat", format="gro")   # treat as GROMACS GRO

# Via read() on an existing object
coor = Coor()
coor.read("structure.dat", format="cif")
```

Accepted `format` values: `"pdb"`, `"cif"`, `"pqr"`, `"gro"`.
`format` defaults to `""`, which means *infer from extension* — the same
behaviour as before.

### Fetching from the PDB archive

```python
coor = Coor(pdb_id="1y0m")    # downloads mmCIF, caches in /tmp/pdb_cpp_cache/
```

The cached file is reused on subsequent calls with the same PDB ID.

---

## Interaction analysis

Interface-oriented tools are grouped under `pdb_cpp.analysis.interaction`:

- `interaction.hbonds()` for hydrogen bonds
- `interaction.salt_bridges()` for ionic contacts / salt bridges
- `interaction.interface_sasa()` for buried surface and interface area

Use `pdb_cpp.analysis.interaction` for interface-oriented analysis. The older
top-level SASA and salt-bridge modules were removed; interface-oriented helpers
live under `pdb_cpp.analysis`.

```python
from pdb_cpp import Coor
from pdb_cpp.analysis import interaction

coor = Coor("tests/input/1A0A.cif")

hb = interaction.hbonds(coor, donor_sel="protein", acceptor_sel="nucleic")
sb = interaction.salt_bridges(coor, cation_sel="protein", anion_sel="nucleic")
iface = interaction.interface_sasa(
    coor,
    receptor_sel="chain C D",
    ligand_sel="chain A B",
)[0]

print(len(hb[0]), len(sb[0]), iface["interface_area"])
```

### Loading through the `rcsb` module

The `pdb_cpp.rcsb` module gives explicit control over which RCSB coordinate
file is downloaded.

```python
from pdb_cpp import rcsb

# Deposited asymmetric unit
asym_unit = rcsb.load("1y0m", structure="asymmetric_unit")

# Biological assembly 1 in mmCIF format
assembly_1 = rcsb.load("5a9z", structure="biological_assembly", assembly_id=1)

# Download without loading, using a custom cache directory
local_path = rcsb.download(
    "5a9z",
    structure="biological_assembly",
    assembly_id=1,
    cache_dir="/tmp/pdb_cpp_rcsb",
)
```

Supported `structure` values are:

- `"asymmetric_unit"`: the deposited entry file, for example `1y0m.cif`
- `"biological_assembly"`: an assembly-specific RCSB file, for example `5a9z-assembly1.cif`

Supported `file_format` values are `"cif"` and `"pdb"`. For biological
assemblies, the RCSB URL pattern differs by format:

- mmCIF: `<pdb_id>-assembly<assembly_id>.cif`
- PDB: `<pdb_id>.pdb<assembly_id>`

If you prefer the older convenience API, `Coor(pdb_id="...")` still works and
now accepts extra RCSB-related keyword arguments:

```python
from pdb_cpp import Coor

coor = Coor(
    pdb_id="5a9z",
    rcsb_structure="biological_assembly",
    assembly_id=1,
    cache_dir="/tmp/pdb_cpp_rcsb",
    force_download=False,
)
```

This path currently reads the downloaded file as returned by RCSB. It does not
construct assemblies locally from `_pdbx_struct_assembly_gen`; instead it relies
on the explicit assembly files published by RCSB.

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
| `conect`          | `dict[int, list[int]]` | CONECT bond records (atom number → bonded atom numbers) |

### `Coor` methods

| Method                              | Returns        | Description                                       |
|-------------------------------------|----------------|---------------------------------------------------|
| `select_atoms(selection)`           | `Coor`         | Select atoms by string query (see §3)             |
| `select_bool_index(bool_array)`     | `Coor`         | Select atoms by boolean mask (one bool per atom)  |
| `get_index_select(selection)`       | `list[int]`    | Atom indices matching a selection string          |
| `read(filename)`                    | `bool`         | Read a file into an existing `Coor` object        |
| `write(filename)`                   | —              | Write to PDB/mmCIF/PQR/GRO (by extension)        |
| `add_Model(model)`                  | —              | Append a `Model` to this `Coor`                   |
| `clear()`                           | —              | Remove all models                                 |
| `get_uniq_chain()`                  | `list`         | Unique chain IDs (char arrays)                    |
| `get_uniq_chain_str()`              | `list[str]`    | Unique chain IDs as Python strings                |
| `get_aa_seq(gap_in_seq=True)`       | `dict`         | Amino acid sequences per chain                    |
| `get_aa_DL_seq()`                   | `dict`         | D/L amino acid sequences (D as lowercase)         |
| `get_aa_na_seq()`                   | `dict`         | Amino acid + nucleic acid sequences               |
| `remove_incomplete_backbone_residues()` | `Coor`     | Remove residues missing backbone atoms            |

### `Model` methods

| Method                              | Returns         | Description                                      |
|-------------------------------------|-----------------|--------------------------------------------------|
| `select_atoms(selection)`           | `list[bool]`    | Boolean mask of atoms matching selection          |
| `sasa(...)`                         | `dict`          | SASA total with optional per-atom/per-residue data |
| `addAtom(...)`                      | `bool`          | Add a single atom (see below)                    |
| `clear()`                           | —               | Remove all atoms from this model                 |
| `get_centroid(indices=None)`        | `ndarray(3,)`   | Centroid of all atoms, or of given indices        |

Note that `Model.select_atoms()` returns a **boolean mask** (not a new
object), unlike `Coor.select_atoms()` which returns a new `Coor`.
`Model.addAtom()` takes positional arguments:
`num, name, resname, resid, chain, x, y, z, occ, beta, altloc, elem, insertres, field, uniq_resid`.

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

### Boolean-mask selection

When you already have a boolean array (e.g. from NumPy operations),
use `select_bool_index` instead of building a string query:

```python
import numpy as np

coor = Coor("tests/input/1y0m.cif")
mask = np.array(coor.beta) < 30.0          # low B-factor atoms
low_b = coor.select_bool_index(mask.tolist())
print(f"{low_b.len} atoms with B < 30")
```

---

## 3.1 SASA and interface SASA

`sasa.sasa()` computes solvent-accessible surface area for each model in a
`Coor` object. It returns one dictionary per model, each with at least
`total`, `polar`, and `apolar`, and optionally `atom_areas` and
`residue_areas`.

### SASA for one chain or selection

```python
from pdb_cpp import Coor
from pdb_cpp.analysis import sasa

coor = Coor("tests/input/1a2k.pdb")
result = sasa.sasa(coor, selection="chain A", by_residue=True)[0]
print(f"Total SASA: {result['total']:.2f} A^2")
print(f"Polar SASA: {result['polar']:.2f} A^2")
print(f"Apolar SASA: {result['apolar']:.2f} A^2")
print(result["residue_areas"][:3])
```

The polar/apolar split is FreeSASA-like: atom areas are partitioned from the
per-atom SASA using a simple element-based classifier, with `N`, `O`, `S`,
`P`, and `Se` treated as polar atoms.

### Protein-protein interface SASA

Use `pdb_cpp.analysis.sasa.buried_surface_area()` to evaluate two partners inside one
complex. This helper computes three SASA values internally:

- receptor alone
- ligand alone
- receptor + ligand together as the complex

The returned interface metrics follow the standard convention:

- `buried_surface = sasa_receptor + sasa_ligand - sasa_complex`
- `interface_area = buried_surface / 2`

```python
from pdb_cpp import Coor
from pdb_cpp.analysis import sasa

coor = Coor("tests/input/1a2k.pdb")

interface = sasa.buried_surface_area(
    coor,
    receptor_sel="chain A",
    ligand_sel="chain B",
    by_residue=True,
)[0]

print(f"Complex SASA   : {interface['complex_sasa']:.2f} A^2")
print(f"Complex polar  : {interface['complex_polar_sasa']:.2f} A^2")
print(f"Complex apolar : {interface['complex_apolar_sasa']:.2f} A^2")
print(f"Receptor SASA  : {interface['receptor_sasa']:.2f} A^2")
print(f"Ligand SASA    : {interface['ligand_sasa']:.2f} A^2")
print(f"Buried surface : {interface['buried_surface']:.2f} A^2")
print(f"Buried polar   : {interface['buried_polar_surface']:.2f} A^2")
print(f"Buried apolar  : {interface['buried_apolar_surface']:.2f} A^2")
print(f"Interface area : {interface['interface_area']:.2f} A^2")
print(f"Interface polar: {interface['interface_polar_area']:.2f} A^2")
print(f"Interface apol.: {interface['interface_apolar_area']:.2f} A^2")

for residue in interface["residue_buried_surface"]:
    print(
        residue["partner"],
        residue["chain"],
        residue["resname"],
        residue["resid"],
        residue["buried_area"],
    )
```

The `residue_buried_surface` list contains one entry per residue from the
receptor and ligand selections, with:

- `isolated_area`: SASA of that residue in the isolated partner
- `complex_area`: SASA of that residue inside the full complex
- `buried_area`: difference between the two

For residue-level SASA, each residue entry also contains `polar_area` and
`apolar_area`; interface residue burial entries additionally expose
`buried_polar_area` and `buried_apolar_area`.

### CONECT records

CONECT records (covalent bonds) are stored as a dictionary mapping atom
numbers to lists of bonded atom numbers. Bond topology is read from and
written to both supported formats:

- **PDB**: `CONECT` lines (hybrid-36 encoded atom serial numbers)
- **mmCIF**: `_struct_conn` loop (any `conn_type_id` is accepted on read;
  `covale` is written on output)

```python
from pdb_cpp import Coor

# Load a structure with covalent bonds (ligands, zinc coordination, etc.)
coor = Coor("tests/input/1u85.pdb")

# Inspect bonds
for atom_num, bonded in list(coor.conect.items())[:3]:
    print(f"Atom {atom_num} bonded to {bonded}")

# Bonds survive PDB and mmCIF round-trips
coor.write("output_with_conect.pdb")   # writes CONECT lines
coor.write("output_with_conect.cif")   # writes _struct_conn loop
```

```{note}
mmCIF files distributed by the PDB typically store only inter-residue bonds
in `_struct_conn` (disulfides, metal coordination, covalent cross-links).
Intra-ligand bonds may be absent unless the file was produced by pdb_cpp.
```

### Programmatic construction

You can build a `Coor` from scratch or append models:

```python
from pdb_cpp import Coor

coor_1 = Coor("tests/input/1y0m.cif")
coor_2 = Coor("tests/input/1y0m.cif")

# Append all models from coor_2 into coor_1
for model in coor_2.models:
    coor_1.add_Model(model)

print(coor_1.model_num)  # doubled

# Clear everything
coor_1.clear()
print(coor_1.model_num)  # 0
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
| `nucleic`  | RNA and DNA residue names                    |
| `rna`      | RNA residue names (`A U G C`)               |
| `dna`      | DNA residue names (`DA DC DG DT`)           |

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

### Low-level C++ alignment: `core.sequence_align`

The `alignment.align_seq()` function is a Python wrapper around the C++
`core.sequence_align()` function, which returns an `Alignment_cpp` object:

```python
from pdb_cpp import core

result = core.sequence_align("ACDEFG", "ACDEG")
print(result.seq1)    # aligned sequence 1 (with gaps)
print(result.seq2)    # aligned sequence 2 (with gaps)
print(result.score)   # alignment score (int)
```

#### `Alignment_cpp` fields

| Field   | Type  | Description                      |
|---------|-------|----------------------------------|
| `seq1`  | `str` | Aligned sequence 1 (with gaps)   |
| `seq2`  | `str` | Aligned sequence 2 (with gaps)   |
| `score` | `int` | Raw alignment score              |

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

The analysis API is also grouped into submodules, so both of these patterns
are supported:

```python
from pdb_cpp import analysis
from pdb_cpp.analysis import dockq, sasa, hbonds
```

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

### Low-level: `core.compute_SS`

The `TMalign.compute_secondary_structure()` wraps the C++ function
`core.compute_SS()`, which returns secondary structure strings directly:

```python
from pdb_cpp import Coor, core

coor = Coor("tests/input/1y0m.cif")
ss = core.compute_SS(coor, gap_in_seq=False)
# Returns list[list[str]]: one list of SS strings per model
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

### Command-line multimer scoring

The installed `pdb_cpp_dockq` executable runs `analysis.dockQ_multimer()`.
This is the same multimer code path used in the benchmark script, including
automatic native-to-model chain-map discovery when `--chain-map` is omitted.

```bash
pdb_cpp_dockq tests/input/1a2k_model.pdb tests/input/1a2k.pdb
```

To provide an explicit mapping, pass native:model pairs:

```bash
pdb_cpp_dockq tests/input/1a2k_model.pdb tests/input/1a2k.pdb --chain-map A:B,B:A,C:C
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

#### `native_contact` parameters

| Parameter                | Default | Description                                        |
|--------------------------|---------|----------------------------------------------------|
| `coor`                   | —       | Model `Coor`                                       |
| `native_coor`            | —       | Native `Coor`                                      |
| `rec_chains`             | —       | Model receptor chain IDs                            |
| `lig_chains`             | —       | Model ligand chain IDs                              |
| `native_rec_chains`      | —       | Native receptor chain IDs                           |
| `native_lig_chains`      | —       | Native ligand chain IDs                             |
| `cutoff`                 | `5.0`   | Contact distance cutoff (Å)                         |
| `residue_id_map`         | `None`  | `dict[int, int]` mapping model residue IDs to a shared ID space |
| `native_residue_id_map`  | `None`  | `dict[int, int]` mapping native residue IDs to the same shared ID space |

The `residue_id_map` parameters are useful when model and native have
different residue numbering — provide dictionaries that map each residue
`uniq_resid` to a common numbering scheme.

#### `interface_rmsd` notes

`interface_rmsd` returns a list of floats with one entry per model. Entries
are `None` when no interface residues are found within the cutoff distance.


### Citation

If you use DockQ scoring, please cite:

- Mirabello, C. & Wallner, B. (2024) *Bioinformatics*. DOI: 10.1093/bioinformatics/btae586

---

## 10. Hydrogen bond detection

`pdb_cpp.analysis.hbonds` detects hydrogen bonds using a purely geometric approach
(Baker & Hubbard, 1984). For structures without explicit hydrogen atoms,
backbone NH and sidechain hydrogen positions are reconstructed
algebraically. The algorithm covers all standard L/D amino acids and RNA/DNA
nucleotides.

### Basic usage

```python
from pdb_cpp import Coor
from pdb_cpp.analysis import hbonds

coor = Coor("tests/input/2rri.cif")

# All protein–protein H-bonds for every model
all_bonds = hbonds.hbonds(coor)

# List of HBond objects for model 0
bonds_model0 = all_bonds[0]
print(f"Model 0: {len(bonds_model0)} hydrogen bonds")

for b in bonds_model0[:3]:
    print(b)
```

### Function signature

```python
hbonds.hbonds(
    coor,
    donor_sel    = "protein",
    acceptor_sel = "protein",
    dist_DA_cutoff = 3.5,   # D···A distance cutoff (Å)
    dist_HA_cutoff = 2.5,   # H···A distance cutoff (Å)
    angle_cutoff   = 90.0,  # D−H···A angle cutoff (°)
) -> list[list[HBond]]
```

Returns a list with one inner list (of `HBond` objects) per model frame.

#### Parameters

| Parameter        | Default     | Description                                           |
|------------------|-------------|-------------------------------------------------------|
| `coor`           | —           | A `Coor` object                                       |
| `donor_sel`      | `"protein"` | Selection string for donor-containing atoms           |
| `acceptor_sel`   | `"protein"` | Selection string for acceptor-containing atoms        |
| `dist_DA_cutoff` | `3.5`       | Maximum heavy-atom donor–acceptor distance (Å)        |
| `dist_HA_cutoff` | `2.5`       | Maximum hydrogen–acceptor distance (Å)                |
| `angle_cutoff`   | `90.0`      | Minimum D−H···A angle (°)                             |

### `HBond` object

Each entry in the returned lists is an `HBond` object with the following
read-only attributes:

| Attribute            | Type    | Description                                        |
|----------------------|---------|----------------------------------------------------|
| `donor_resid`        | `int`   | Donor residue sequence number                      |
| `donor_resname`      | `str`   | Donor residue name                                 |
| `donor_chain`        | `str`   | Donor chain identifier                             |
| `donor_heavy_name`   | `str`   | Donor heavy-atom name (e.g. `N`, `OG`)             |
| `donor_h_name`       | `str`   | Donor hydrogen atom name (e.g. `H`, `HG`)         |
| `donor_heavy_xyz`    | `tuple` | Donor heavy-atom coordinates (x, y, z)             |
| `donor_h_xyz`        | `tuple` | Donor hydrogen coordinates (x, y, z)               |
| `acceptor_resid`     | `int`   | Acceptor residue sequence number                   |
| `acceptor_resname`   | `str`   | Acceptor residue name                              |
| `acceptor_chain`     | `str`   | Acceptor chain identifier                          |
| `acceptor_name`      | `str`   | Acceptor atom name (e.g. `O`, `OD1`)              |
| `acceptor_xyz`       | `tuple` | Acceptor atom coordinates (x, y, z)                |
| `dist_DA`            | `float` | D···A distance (Å)                                 |
| `dist_HA`            | `float` | H···A distance (Å)                                 |
| `angle_DHA`          | `float` | D−H···A angle (°)                                  |

### Filtering and analysis examples

```python
from pdb_cpp import Coor
from pdb_cpp.analysis import hbonds

coor = Coor("tests/input/2rri.cif")
bonds = hbonds.hbonds(coor)[0]

# Backbone–backbone H-bonds only
bb_names = {"N", "O"}
bb_bonds = [b for b in bonds
            if b.donor_heavy_name in bb_names
            and b.acceptor_name in bb_names]

# Inter-chain H-bonds (cross-chain contacts)
inter = [b for b in bonds if b.donor_chain != b.acceptor_chain]
print(f"{len(inter)} inter-chain H-bonds")

# Sort by D–A distance
bonds.sort(key=lambda b: b.dist_DA)
for b in bonds[:5]:
    print(f"{b.donor_chain}{b.donor_resid}{b.donor_resname:>4}"
          f" {b.donor_heavy_name} → "
          f"{b.acceptor_chain}{b.acceptor_resid}{b.acceptor_resname:>4}"
          f" {b.acceptor_name}  d(DA)={b.dist_DA:.2f} Å"
          f"  angle={b.angle_DHA:.1f}°")
```

### Cross-selection H-bonds (e.g. protein–nucleic)

The `donor_sel` and `acceptor_sel` parameters accept any valid selection
string, enabling inter-molecular or cross-type analysis:

```python
# Protein side-chains donating to RNA acceptors
rna_contacts = hbonds.hbonds(
    coor,
    donor_sel    = "protein",
    acceptor_sel = "nucleic",
)

# H-bonds between two specific chains
ab_bonds = hbonds.hbonds(coor, donor_sel="chain A", acceptor_sel="chain B")
```

### Baker & Hubbard criteria

The default thresholds reproduce the original Baker & Hubbard (1984)
geometric criteria:

$$d(\text{D}{\cdots}\text{A}) \leq 3.5~\text{Å}, \quad
  d(\text{H}{\cdots}\text{A}) \leq 2.5~\text{Å}, \quad
  \angle(\text{D{-}H}{\cdots}\text{A}) \geq 90°$$

A stricter cut-off of 120° (`angle_cutoff=120`) is closer to values used
by biotite and is recommended when comparing with other tools.

---

## 11. Geometry utilities

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

## 12. Data cleaning utilities

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

## 13. Module summary

| Module              | Key functions / classes                                   |
|---------------------|-----------------------------------------------------------|
| `pdb_cpp`           | `Coor`, `Model` — main data objects                           |
| `pdb_cpp.core`      | `get_common_atoms`, `coor_align`, `align_seq_based`, `tmalign_ca`, `distance_matrix`, `compute_SS`, `sequence_align`, `Alignment_cpp`, `hy36encode`, `hy36decode` |
| `pdb_cpp.alignment` | `align_seq`, `print_align_seq`, `align_chain_permutation` |
| `pdb_cpp.analysis`  | `rmsd`, `interface_rmsd`, `native_contact`, `dockQ`       |
| `pdb_cpp.analysis.hbonds` | `hbonds` — Baker & Hubbard H-bond detection         |
| `pdb_cpp.TMalign`   | `compute_secondary_structure`                             |
| `pdb_cpp.geom`      | `distance_matrix`                                         |
| `pdb_cpp.select`    | `remove_incomplete_backbone_residues`                     |
| `pdb_cpp.sequence`  | `get_aa_seq`, `get_aa_DL_seq`                             |
| `pdb_cpp.data`      | `BLOSUM62`, `get_blosum62`, `load_blosum`, `AA_DICT`, `NA_DICT`, `AA_NA_DICT` |

---

## 14. Data subpackage (`pdb_cpp.data`)

The `pdb_cpp.data` module provides residue dictionaries and the BLOSUM62
substitution matrix used internally by the alignment routines.

### BLOSUM62 matrix

```python
from pdb_cpp.data.blosum import BLOSUM62

# Substitution score for two amino acids
print(BLOSUM62[("A", "A")])   # 4
print(BLOSUM62[("A", "W")])   # -3

# Load a custom matrix from file
from pdb_cpp.data.blosum import load_blosum
custom = load_blosum("/path/to/matrix.txt")
```

`BLOSUM62` is a lazy-loaded dictionary mapping `(aa1, aa2)` tuples to
integer substitution scores.

### Residue dictionaries

```python
from pdb_cpp.data.res_dict import AA_DICT, AA_DICT_L, AA_DICT_D, NA_DICT

# 3-letter to 1-letter amino acid code
print(AA_DICT["ALA"])   # "A"
print(AA_DICT["GLY"])   # "G"

# D-amino acids (e.g. DAL -> "A")
print(AA_DICT_D["DAL"])  # "A"

# Nucleic acids
print(NA_DICT["DA"])     # "A"

# Combined amino acid + nucleic acid dictionary
from pdb_cpp.data.res_dict import AA_NA_DICT
```

| Dictionary   | Contents                                    |
|-------------|---------------------------------------------|
| `AA_DICT_L` | Standard L-amino acids (31 entries incl. protonation variants) |
| `AA_DICT_D` | D-amino acids (22 entries)                   |
| `AA_DICT`   | Combined L + D amino acids                   |
| `NA_DICT`   | DNA nucleotides (DA, DT, DC, DG)             |
| `AA_NA_DICT`| Combined amino acids + nucleotides           |
