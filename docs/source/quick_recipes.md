# Quick Recipes

Short, copy-paste snippets for common workflows.

## 1) Load a structure from local file or PDB ID

```python
from pdb_cpp import Coor

coor_local = Coor("tests/input/1y0m.cif")
coor_remote = Coor(pdb_id="1y0m")
```

## 2) Extract an interface selection between two chains

```python
from pdb_cpp import Coor

coor = Coor("tests/input/1rxz.pdb")

# Atoms in chain A within 10 Å of chain B
interface_a = coor.select_atoms("chain A and within 10.0 of chain B")
interface_a.write("interface_A_vs_B.pdb")
```

## 3) Select a receptor-ligand complex subset

```python
from pdb_cpp import Coor

coor = Coor("tests/input/1rxz.pdb")

# Keep only chains A and B (example receptor/ligand pair)
complex_ab = coor.select_atoms("chain A B")
complex_ab.write("complex_AB.pdb")
```

## 4) Sequence alignment for two chains

```python
from pdb_cpp import Coor, alignment

coor_1 = Coor("tests/input/1u85.pdb")
coor_2 = Coor("tests/input/1ubd.pdb")

seq_1 = coor_1.get_aa_seq()["A"]
seq_2 = coor_2.get_aa_seq()["C"]

aln_1, aln_2, score = alignment.align_seq(seq_1, seq_2)
print(score)
alignment.print_align_seq(aln_1, aln_2)
```

## 5) Align coordinates and compute RMSD

```python
from pdb_cpp import Coor, core, analysis

coor_1 = Coor("tests/input/1u85.pdb")
coor_2 = Coor("tests/input/1ubd.pdb")

idx_1, idx_2 = core.get_common_atoms(coor_1, coor_2, chain_1=["A"], chain_2=["C"])
core.coor_align(coor_1, coor_2, idx_1, idx_2, frame_ref=0)

rmsd_values = analysis.rmsd(coor_1, coor_2, index_list=[idx_1, idx_2])
print(rmsd_values[0])
```

## 6) TM-score for specific chain pair

```python
from pdb_cpp import Coor
from pdb_cpp.core import tmalign_ca

coor_1 = Coor("tests/input/1y0m.cif")
coor_2 = Coor("tests/input/1ubd.pdb")

tm = tmalign_ca(coor_1, coor_2, chain_1=["A"], chain_2=["C"], mm=1)
print(tm.L_ali, tm.rmsd, tm.TM1, tm.TM2)
```

If you use this USalign/TM-align functionality, please cite:

- Chengxin Zhang, Morgan Shine, Anna Marie Pyle, Yang Zhang (2022) Nat Methods. 19(9), 1109-1115.
- Chengxin Zhang, Anna Marie Pyle (2022) iScience. 25(10), 105218.

## 7) DockQ with automatic chain-role inference

```python
from pdb_cpp import Coor, analysis

model = Coor("tests/input/1rxz_colabfold_model_1.pdb")
native = Coor("tests/input/1rxz.pdb")

scores = analysis.dockQ(model, native)
print(scores["DockQ"][0])
print(scores["Fnat"][0], scores["Fnonnat"][0])
print(scores["LRMS"][0], scores["iRMS"][0], scores["rRMS"][0])
```

## 8) DockQ with explicit receptor/ligand chains

```python
from pdb_cpp import Coor, analysis

model = Coor("tests/input/1rxz_colabfold_model_1.pdb")
native = Coor("tests/input/1rxz.pdb")

scores = analysis.dockQ(
    model,
    native,
    rec_chains=["B"],
    lig_chains=["C"],
    native_rec_chains=["A"],
    native_lig_chains=["B"],
)
print(scores["DockQ"][0])
```

If you use DockQ scoring, please cite:

- DockQ, DOI: 10.1093/bioinformatics/btae586

## 9) Secondary structure per model/chain

```python
from pdb_cpp import Coor, TMalign

coor = Coor("tests/input/1y0m.cif")
ss_list = TMalign.compute_secondary_structure(coor)

print(ss_list[0]["A"])
```

## 10) Distance matrix on C-alpha atoms

```python
from pdb_cpp import Coor, geom

coor = Coor("tests/input/1y0m.cif")
ca = coor.select_atoms("name CA")
dmat = geom.distance_matrix(ca, ca)
print(dmat.shape)
```
