# Basic Usage

A quick-start walk-through. For full details on every feature see the
[Functionality Guide](functionality.md); for copy-paste snippets see
[Quick Recipes](quick_recipes.md).

## Import and load structures

```python
from pdb_cpp import Coor

# Load from local file (PDB, mmCIF, PQR, or GRO — auto-detected by extension)
coor = Coor("tests/input/1y0m.cif")

# Load from PDB ID (downloaded as mmCIF and cached locally)
coor_from_id = Coor(pdb_id="1y0m")

# Load explicitly from the RCSB helper
from pdb_cpp import rcsb

asym_unit = rcsb.load("1y0m", structure="asymmetric_unit")
bio_assembly = rcsb.load("5a9z", structure="biological_assembly", assembly_id=1)
```

Use `Coor(pdb_id=...)` when you only need the deposited mmCIF entry. Use
`pdb_cpp.rcsb` when you want explicit control over which RCSB file is fetched,
including biological assemblies, cache location, or forced re-download.

## Inspect a structure

```python
print(coor.model_num)            # number of models
print(coor.len)                  # number of atoms in active model
print(coor.get_uniq_chain())     # list of chain IDs

# Coordinates, atom names, residue names
print(coor.xyz[:5])              # shape (N, 3)
print(coor.name_str[:5])         # ['N', 'CA', 'C', 'O', 'CB']
print(coor.resname_str[:5])      # ['THR', 'THR', ...]
print(coor.chain_str[:5])        # ['A', 'A', ...]
print(coor.resid[:5])            # residue sequence numbers
print(coor.beta[:5])             # B-factors
print(coor.occ[:5])              # occupancies
```

See the [Functionality Guide](functionality.md) for the complete list of
`Coor` and `Model` properties.

## Atom selections

```python
chain_a = coor.select_atoms("chain A")
backbone = coor.select_atoms("protein and backbone")
interface = coor.select_atoms("chain A and within 5.0 of chain B")
indices = coor.get_index_select("name CA and chain A")
```

See the [Functionality Guide](functionality.md) for the full selection syntax,
including keywords, operators, and spatial queries.

## Write output

```python
chain_a.write("chain_A.pdb")
chain_a.write("chain_A.cif")
```

## Sequence alignment and structural superposition

```python
from pdb_cpp import alignment, core, analysis

coor_1 = Coor("tests/input/1u85.pdb")
coor_2 = Coor("tests/input/1ubd.pdb")

# One-step sequence-based alignment + RMSD
rmsd_list, _, _ = core.align_seq_based(coor_1, coor_2, chain_1=["A"], chain_2=["C"])
print(f"RMSD: {rmsd_list[0]:.3f} Å")
```

## TM-score

```python
from pdb_cpp.core import tmalign_ca

tm = tmalign_ca(coor_1, coor_2, chain_1=["A"], chain_2=["C"], mm=1)
print(f"TM1={tm.TM1:.4f}, TM2={tm.TM2:.4f}, RMSD={tm.rmsd:.3f}")
```

## DockQ scoring

```python
model = Coor("tests/input/1rxz_colabfold_model_1.pdb")
native = Coor("tests/input/1rxz.pdb")

scores = analysis.dockQ(model, native)
print(f"DockQ: {scores['DockQ'][0]:.3f}")
```

For multimer scoring from the command line, use the installed
`pdb_cpp_dockq` executable. It calls `analysis.dockQ_multimer()` and uses
the same automatic native-to-model chain mapping as the benchmark script:

```bash
pdb_cpp_dockq tests/input/1a2k_model.pdb tests/input/1a2k.pdb
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

