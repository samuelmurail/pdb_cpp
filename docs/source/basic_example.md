# Basic Usage

## Import and load structures

```python
from pdb_cpp import Coor

# Load from local file
coor = Coor("tests/input/1y0m.cif")

# Load from PDB ID (downloaded as mmCIF and cached locally)
coor_from_id = Coor(pdb_id="1y0m")
```

`Coor` can read `.pdb` and `.cif` files and can contain one or more models.

## Read basic properties

```python
print(coor.model_num)        # number of models
print(coor.len)              # number of atoms in active model
print(coor.get_uniq_chain()) # list of chain IDs
print(coor.get_aa_seq())     # dict: chain -> sequence
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

seq_1 = coor_1.get_aa_seq()["A"]
seq_2 = coor_2.get_aa_seq()["C"]

aln_1, aln_2, score = alignment.align_seq(seq_1, seq_2)
alignment.print_align_seq(aln_1, aln_2)

idx_1, idx_2 = core.get_common_atoms(coor_1, coor_2, chain_1=["A"], chain_2=["C"])
core.coor_align(coor_1, coor_2, idx_1, idx_2, frame_ref=0)

rmsd = analysis.rmsd(coor_1, coor_2, index_list=[idx_1, idx_2])
print(rmsd[0])
```

## TM-score / TM-align

```python
from pdb_cpp.core import tmalign_ca

coor_1 = Coor("tests/input/1y0m.cif")
coor_2 = Coor("tests/input/1ubd.pdb")

tm = tmalign_ca(coor_1, coor_2, chain_1=["A"], chain_2=["C"], mm=1)

print(tm.L_ali, tm.rmsd)
print(tm.TM1, tm.TM2)
```

## DockQ scoring

```python
from pdb_cpp import analysis

model = Coor("tests/input/1rxz_colabfold_model_1.pdb")
native = Coor("tests/input/1rxz.pdb")

scores = analysis.dockQ(model, native)
print(scores["DockQ"][0])
print(scores["Fnat"][0], scores["Fnonnat"][0])
print(scores["LRMS"][0], scores["iRMS"][0], scores["rRMS"][0])
```

## Secondary structure

```python
from pdb_cpp import TMalign

ss = TMalign.compute_secondary_structure(coor)
print(ss[0]["A"])
```

