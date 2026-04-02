# Quick Recipes

Short, copy-paste snippets for common workflows.

---

## 1) Load a structure from local file or PDB ID

```python
from pdb_cpp import Coor

coor_local = Coor("tests/input/1y0m.cif")
coor_remote = Coor(pdb_id="1y0m")
```

## 2) Inspect atom properties

```python
from pdb_cpp import Coor

coor = Coor("tests/input/1y0m.cif")

# Coordinates as a NumPy array
print(coor.xyz.shape)             # (N, 3)

# Atom names, chain IDs, residue names as strings
print(coor.name_str[:5])          # ['N', 'CA', 'C', 'O', 'CB']
print(coor.chain_str[:5])         # ['A', 'A', 'A', 'A', 'A']
print(coor.resname_str[:5])       # ['THR', 'THR', ...]

# B-factors and occupancies
print(coor.beta[:5])
print(coor.occ[:5])
```

## 3) Extract an interface selection between two chains

```python
from pdb_cpp import Coor

coor = Coor("tests/input/1rxz.pdb")

# Atoms in chain A within 10 Å of chain B
interface_a = coor.select_atoms("chain A and within 10.0 of chain B")
interface_a.write("interface_A_vs_B.pdb")
```

## 4) Select a receptor-ligand complex subset

```python
from pdb_cpp import Coor

coor = Coor("tests/input/1rxz.pdb")

# Keep only chains A and B (example receptor/ligand pair)
complex_ab = coor.select_atoms("chain A B")
complex_ab.write("complex_AB.pdb")
```

## 5) Clean up: remove incomplete backbone residues

```python
from pdb_cpp import Coor, select

coor = Coor("tests/input/1y0m.cif")
clean = select.remove_incomplete_backbone_residues(coor)
print(f"Before: {coor.len}, After: {clean.len}")
```

## 6) Sequence alignment for two chains

```python
from pdb_cpp import Coor, alignment

coor_1 = Coor("tests/input/1u85.pdb")
coor_2 = Coor("tests/input/1ubd.pdb")

seq_1 = coor_1.get_aa_seq()["A"]
seq_2 = coor_2.get_aa_seq()["C"]

aln_1, aln_2, score = alignment.align_seq(seq_1, seq_2)
print(f"Score: {score}")
alignment.print_align_seq(aln_1, aln_2)
```

## 7) Align coordinates and compute RMSD

```python
from pdb_cpp import Coor, core, analysis

coor_1 = Coor("tests/input/1u85.pdb")
coor_2 = Coor("tests/input/1ubd.pdb")

idx_1, idx_2 = core.get_common_atoms(coor_1, coor_2, chain_1=["A"], chain_2=["C"])
core.coor_align(coor_1, coor_2, idx_1, idx_2, frame_ref=0)

rmsd_values = analysis.rmsd(coor_1, coor_2, index_list=[idx_1, idx_2])
print(f"RMSD: {rmsd_values[0]:.3f} Å")
```

## 8) One-step sequence-based alignment

```python
from pdb_cpp import Coor, core

coor_1 = Coor("tests/input/1u85.pdb")
coor_2 = Coor("tests/input/1ubd.pdb")

rmsd_list, align_idx_1, align_idx_2 = core.align_seq_based(
    coor_1, coor_2, chain_1=["A"], chain_2=["C"]
)
print(f"RMSD: {rmsd_list[0]:.3f} Å")
```

## 9) Chain-permutation alignment (unknown chain mapping)

```python
from pdb_cpp import Coor, alignment

coor_1 = Coor("tests/input/1u85.pdb")
coor_2 = Coor("tests/input/1ubd.pdb")

rmsds, mappings = alignment.align_chain_permutation(coor_1, coor_2)
print(f"Best RMSD: {rmsds[0]:.3f} Å")
```

## 10) TM-score for specific chain pair

```python
from pdb_cpp import Coor
from pdb_cpp.core import tmalign_ca

coor_1 = Coor("tests/input/1y0m.cif")
coor_2 = Coor("tests/input/1ubd.pdb")

tm = tmalign_ca(coor_1, coor_2, chain_1=["A"], chain_2=["C"], mm=1)
print(f"L_ali={tm.L_ali}, RMSD={tm.rmsd:.3f}, TM1={tm.TM1:.4f}, TM2={tm.TM2:.4f}")
```

If you use this USalign/TM-align functionality, please cite:

- Chengxin Zhang, Morgan Shine, Anna Marie Pyle, Yang Zhang (2022) Nat Methods. 19(9), 1109-1115.
- Chengxin Zhang, Anna Marie Pyle (2022) iScience. 25(10), 105218.

## 11) DockQ with automatic chain-role inference

```python
from pdb_cpp import Coor, analysis

model = Coor("tests/input/1rxz_colabfold_model_1.pdb")
native = Coor("tests/input/1rxz.pdb")

scores = analysis.dockQ(model, native)
print(f"DockQ: {scores['DockQ'][0]:.3f}")
print(f"Fnat: {scores['Fnat'][0]:.3f}, Fnonnat: {scores['Fnonnat'][0]:.3f}")
print(f"LRMS: {scores['LRMS'][0]:.3f}, iRMS: {scores['iRMS'][0]:.3f}")
```

## 12) DockQ with explicit receptor/ligand chains

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
print(f"DockQ: {scores['DockQ'][0]:.3f}")
```

If you use DockQ scoring, please cite:

- DockQ, DOI: 10.1093/bioinformatics/btae586

## 13) Secondary structure per model/chain

```python
from pdb_cpp import Coor, TMalign

coor = Coor("tests/input/1y0m.cif")
ss_list = TMalign.compute_secondary_structure(coor)

for chain_id, ss_string in ss_list[0].items():
    print(f"Chain {chain_id}: {ss_string}")
```

## 14) Distance matrix on C-alpha atoms

```python
from pdb_cpp import Coor, geom

coor = Coor("tests/input/1y0m.cif")
ca = coor.select_atoms("name CA")
dmat = geom.distance_matrix(ca, ca)
print(f"Shape: {dmat.shape}")
```

## 15) D/L amino acid and nucleic acid sequences

```python
from pdb_cpp import Coor

coor = Coor("tests/input/1y0m.cif")

# Standard amino acid sequences
seqs = coor.get_aa_seq()
print(seqs)

# D-residues as lowercase
dl_seqs = coor.get_aa_DL_seq()
print(dl_seqs)

# Amino acid + nucleic acid sequences
all_seqs = coor.get_aa_na_seq()
print(all_seqs)
```

## 16) Hybrid-36 encoding/decoding

```python
from pdb_cpp.core import hy36encode, hy36decode

# Encode large atom numbers for PDB format
encoded = hy36encode(5, 100000)   # "A0000"
decoded = hy36decode(5, "A0000")  # 100000
print(encoded, decoded)
```

## 17) Multi-model: iterate and compute per-model RMSD

```python
from pdb_cpp import Coor, analysis

coor = Coor("tests/input/2rri.pdb")  # NMR ensemble
ref = Coor("tests/input/2rri.pdb")

# RMSD of each model against model 0
rmsd_values = analysis.rmsd(coor, ref, selection="name CA", frame_ref=0)
for i, r in enumerate(rmsd_values):
    print(f"Model {i}: RMSD = {r:.3f} Å")
```

## 18) Bond topology (CONECT / _struct_conn)

```python
from pdb_cpp import Coor

# Both PDB (CONECT lines) and mmCIF (_struct_conn) bond records are supported
coor = Coor("tests/input/1u85.pdb")

# Inspect covalent bonds (dict: atom_serial -> list[bonded_serials])
print(f"Atoms with bonds: {len(coor.conect)}")
for atom_num, bonded in list(coor.conect.items())[:5]:
    print(f"  Atom {atom_num} -> {bonded}")

# Convert PDB to mmCIF, preserving bond topology via _struct_conn
coor.write("output_1u85.cif")

# Reload from mmCIF and confirm bonds survived
coor2 = Coor("output_1u85.cif")
print(f"Bonds after mmCIF round-trip: {len(coor2.conect)}")

# Bonds are also preserved after atom selection and renumbering
ligand = coor.select_atoms("not protein")
ligand.write("ligand_only.pdb")
```

## 19) Hydrogen bond detection

```python
from pdb_cpp import Coor
from pdb_cpp import hbond

coor = Coor("tests/input/2rri.cif")

# All protein–protein H-bonds (one list per model)
all_bonds = hbond.hbonds(coor)
print(f"Model 0: {len(all_bonds[0])} H-bonds")

# Inspect an individual bond
b = all_bonds[0][0]
print(f"{b.donor_chain}{b.donor_resid} {b.donor_heavy_name}"
      f" -> {b.acceptor_chain}{b.acceptor_resid} {b.acceptor_name}"
      f"  d(DA)={b.dist_DA:.2f} Å  angle={b.angle_DHA:.1f}°")

# Filter inter-chain H-bonds
inter = [b for b in all_bonds[0] if b.donor_chain != b.acceptor_chain]
print(f"{len(inter)} inter-chain H-bonds in model 0")

# Protein donors → nucleic-acid acceptors (protein–RNA interface)
rna_bonds = hbond.hbonds(coor, donor_sel="protein", acceptor_sel="nucleic")
```
