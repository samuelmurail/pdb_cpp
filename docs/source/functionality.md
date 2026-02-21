# Functionality Guide

This guide summarizes the main `pdb_cpp` workflows with links to the corresponding modules.

## 1) Read and write structures

Use `Coor` to load and write coordinates:

```python
from pdb_cpp import Coor

coor = Coor("tests/input/1y0m.cif")
coor.write("copy_1y0m.pdb")
```

`Coor` accepts:

- local `.pdb` and `.cif` files
- `pdb_id="XXXX"` to fetch mmCIF from the PDB archive

## 2) Atom selection and complex/interface selection

Selections are done with `Coor.select_atoms(selection, frame=...)`.

```python
# protein atoms from chain A
prot_a = coor.select_atoms("protein and chain A")

# residues near partner chain (interface-like selection)
interface_a = coor.select_atoms("chain A and within 10.0 of chain B")

# explicit residue and atom filters
subset = coor.select_atoms("name CA N C O and residue >= 10 and residue <= 50")
```

Selection language supports:

- fields: `name`, `altloc`, `resname`, `chain`, `resid`, `residue`, `x`, `y`, `z`, `beta`, `occ`
- keywords: `protein`, `backbone`, `noh`
- logic: `and`, `or`, `not`
- spatial: `within`
- comparisons: `==`, `!=`, `<`, `<=`, `>`, `>=`

## 3) Sequence extraction and sequence alignment

```python
from pdb_cpp import alignment

seqs = coor.get_aa_seq()
chain_a = seqs["A"]

aln_1, aln_2, score = alignment.align_seq(chain_a, chain_a)
alignment.print_align_seq(aln_1, aln_2)
```

Related modules:

- `pdb_cpp.sequence`
- `pdb_cpp.alignment`

## 4) Structural alignment (RMSD-based)

```python
from pdb_cpp import Coor, core, analysis

coor_1 = Coor("tests/input/1u85.pdb")
coor_2 = Coor("tests/input/1ubd.pdb")

idx_1, idx_2 = core.get_common_atoms(coor_1, coor_2, chain_1=["A"], chain_2=["C"])
core.coor_align(coor_1, coor_2, idx_1, idx_2, frame_ref=0)

rmsd_values = analysis.rmsd(coor_1, coor_2, index_list=[idx_1, idx_2])
```

If chain correspondence is unknown in multi-chain systems, use:

```python
from pdb_cpp import alignment

rmsds, mappings = alignment.align_chain_permutation(coor_1, coor_2)
```

## 5) TM-align and TM-score

```python
from pdb_cpp.core import tmalign_ca

tm = tmalign_ca(coor_1, coor_2, chain_1=["A"], chain_2=["C"], mm=1)
print(tm.L_ali, tm.rmsd, tm.TM1, tm.TM2)
```

Outputs include aligned length, RMSD, and both TM-score normalizations.

If you use the USalign/TM-align functionality in `pdb_cpp`, please cite:

- Chengxin Zhang, Morgan Shine, Anna Marie Pyle, Yang Zhang (2022) Nat Methods. 19(9), 1109-1115.
- Chengxin Zhang, Anna Marie Pyle (2022) iScience. 25(10), 105218.

## 6) DockQ for docking quality

```python
from pdb_cpp import analysis

dockq = analysis.dockQ(coor_1, coor_2)
print(dockq["DockQ"][0])
print(dockq["Fnat"][0], dockq["Fnonnat"][0], dockq["LRMS"][0], dockq["iRMS"][0])
```

`analysis.dockQ` computes:

- `Fnat` / `Fnonnat`
- ligand RMSD (`LRMS`)
- interface RMSD (`iRMS`)
- receptor RMSD (`rRMS`)
- combined DockQ score (`DockQ`)

If you use DockQ scoring in `pdb_cpp`, please cite:

- DockQ, DOI: 10.1093/bioinformatics/btae586

## 7) Secondary structure assignment

```python
from pdb_cpp import TMalign

ss_list = TMalign.compute_secondary_structure(coor)
print(ss_list[0]["A"])
```

## 8) Geometry helper functions

```python
from pdb_cpp import geom

ca = coor.select_atoms("name CA")
dmat = geom.distance_matrix(ca, ca)
```

## 9) Benchmarks and performance scripts

Benchmark scripts are available in `benchmark/` for:

- DockQ comparison
- I/O speed comparison
- common operation comparisons against other libraries
