#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_DIR="${VENV_DIR:-$ROOT_DIR/.venv-vg}"
LOOPS="${LOOPS:-2000}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

cd "$ROOT_DIR"

if [[ ! -d "$VENV_DIR" ]]; then
  "$PYTHON_BIN" -m venv "$VENV_DIR"
fi

# shellcheck disable=SC1090
source "$VENV_DIR/bin/activate"

python -m pip install -U pip setuptools wheel
python -m pip install pybind11==2.13 numpy
CXXFLAGS='-O0 -g -fno-omit-frame-pointer' python -m pip install -e . --no-build-isolation --force-reinstall

PYTHONMALLOC=malloc valgrind \
  --tool=memcheck \
  --leak-check=full \
  --show-leak-kinds=definite,indirect \
  --errors-for-leak-kinds=definite,indirect \
  --track-origins=yes \
  --num-callers=40 \
  python - <<PY
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.getcwd(), "src", "pdb_cpp"))
import core

Coor = core.Coor
probe = os.path.join(tempfile.gettempdir(), "vg_short_probe_core_direct.pdb")
with open(probe, "w", encoding="utf-8") as handle:
    handle.write("\n".join([
        "ATOM      1  CA  ALA A   1      10.000  11.000  12.000  1.00 36.12",
        "ATOM      2  P   DA  B   2      13.000  14.000  15.000  1.00 39.41",
        "HETATM    3  ZN  ZN  C   3      16.000  17.000  18.000  1.00 41.03",
        "HETATM    4  O1  LIG D   4      19.000  20.000  21.000  1.00 42.34",
        "END",
    ]) + "\n")

sel = "(protein and name CA) or (dna and name P) or ions or (not protein and not dna and noh)"
for _ in range(int("$LOOPS")):
    coor = Coor(probe)
    selected = coor.select_atoms(sel)
    _ = selected.size()
    _ = coor.get_Models(0).get_beta()

print("done")
PY
