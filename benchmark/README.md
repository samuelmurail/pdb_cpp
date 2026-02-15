# Benchmark scripts

This folder contains scripts to benchmark `pdb_cpp` against other libraries (`pdb_numpy`, `biopython`, `biotite`) and generate plots.

## Prerequisites

Run from the repository root:

```bash
cd /home/murail/Documents/Code/pdb_cpp
```

Use your benchmark environment (for example your `biopython` conda env), with required packages installed.

## 1) Common operations benchmark (recommended)

Measures:
- `read`
- `write`
- `select_within10_chainA`
- `get_aa_seq`
- `rmsd_ca_shift`
- `dihedral_ca`
- `align_seq_chainA`
- `align_ca_self`

### Run on default files

```bash
python benchmark/compare_common_speed.py --runs 3 --warmup 1
```

Output CSV:
- `benchmark/common_speed_comparison.csv`

### Run on a specific structure (example: 9X0F)

```bash
python benchmark/compare_common_speed.py \
  --runs 5 --warmup 1 \
  --files tests/input/9X0F.cif \
  --csv benchmark/common_speed_9X0F.csv
```

## 2) Build histogram / grouped chart

Use the plotting script on any benchmark CSV:

```bash
python benchmark/plot_io_histogram.py \
  --input benchmark/common_speed_comparison.csv \
  --output benchmark/common_speed_histogram.png
```

For 9X0F:

```bash
python benchmark/plot_io_histogram.py \
  --input benchmark/common_speed_9X0F.csv \
  --output benchmark/common_speed_9X0F.png
```

Notes:
- Y axis is logarithmic by default.
- Use `--linear` for linear scale.
- Legend is placed outside the plot (right side).

## 3) IO-only benchmark (read/write only)

```bash
python benchmark/compare_io_speed.py --runs 3 --warmup 1
```

Output CSV:
- `benchmark/io_speed_comparison.csv`

Plot:

```bash
python benchmark/plot_io_histogram.py \
  --input benchmark/io_speed_comparison.csv \
  --output benchmark/io_speed_histogram.png
```

## 4) DockQ benchmark

```bash
python benchmark/compare_dockq_speed.py --runs 3 --warmup 1 --mode end-to-end
```

Output CSV:
- `benchmark/dockq_vs_pdb_cpp.csv`

Optional markdown table:

```bash
python benchmark/csv_to_markdown.py \
  --input benchmark/dockq_vs_pdb_cpp.csv \
  --output benchmark/dockq_vs_pdb_cpp.md
```

## Quick 9X0F workflow

```bash
python benchmark/compare_common_speed.py \
  --runs 5 --warmup 1 \
  --files tests/input/9X0F.cif \
  --csv benchmark/common_speed_9X0F.csv

python benchmark/plot_io_histogram.py \
  --input benchmark/common_speed_9X0F.csv \
  --output benchmark/common_speed_9X0F.png
```
