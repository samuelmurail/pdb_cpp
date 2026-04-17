# Contributing

## Development setup

```bash
git clone https://github.com/samuelmurail/pdb_cpp.git
cd pdb_cpp
pip install -r requirements.txt
pip install -e . --no-build-isolation
```

## Project layout

```
src/pdb_cpp/
├── __init__.py          # Package entry point (Coor, Model)
├── _pyprops.py          # Python property patches on Coor/Model
├── alignment.py         # Sequence alignment wrappers
├── analysis/            # High-level analysis package
│   ├── __init__.py      # Flat compatibility exports + grouped submodules
│   ├── dockq.py         # RMSD, DockQ, interface metrics
│   ├── sasa.py          # SASA namespace bridge
│   └── hbonds.py        # H-bond namespace bridge
├── geom.py              # Distance matrix wrapper
├── select.py            # Backbone cleaning utility
├── sequence.py          # Sequence extraction wrappers
├── TMalign.py           # Secondary structure wrapper
├── _core/               # C++ extension source
│   ├── pybind.cpp       # pybind11 bindings (pdb_cpp.core)
│   ├── Coor.cpp/h       # Multi-model coordinate container
│   ├── Model.cpp/h      # Single-model atom storage
│   ├── align.cpp/h      # Structural alignment algorithms
│   ├── select.cpp/h     # Selection language parser
│   ├── seq_align.cpp/h  # Needleman-Wunsch alignment
│   ├── sequence.cpp/h   # Sequence extraction
│   ├── geom.h           # Geometry (Kabsch, distance matrix)
│   ├── TMalign_wrapper.cpp  # USalign integration
│   ├── format/          # PDB/mmCIF parsers and writers
│   └── usalign/         # Vendored USalign headers
└── data/                # BLOSUM62 matrix, residue dictionaries
```

## Adding a C++ feature

1. Add `.cpp`/`.h` files in `src/pdb_cpp/_core/`.
2. Register new `.cpp` files in `setup.py`'s `ext_modules` source list.
3. Expose Python bindings in `src/pdb_cpp/_core/pybind.cpp`.
4. Add a Python wrapper in the appropriate `src/pdb_cpp/*.py` module.
5. Add tests in `tests/`.
6. Rebuild: `pip install -e . --no-build-isolation`
7. Run tests: `pytest`

## Running tests

```bash
pytest                     # full suite
pytest tests/test_mmcif.py # single file
pytest -v                  # verbose output
pytest -k "tmalign"        # filter by name
```

## Code style

- Python: follow PEP 8; use type hints in new code.
- C++: use C++17; prefer `std::runtime_error` over `assert()` for
  error conditions that can occur with user input.
- Docstrings: use NumPy-style docstrings for all public functions.

## Building documentation

```bash
cd docs
pip install -r requirements.txt
make html
```

The output is in `docs/build/html/`.

## Memory safety checks

Sanitizer and Valgrind scripts are in `scripts/`:

```bash
bash scripts/asan_core_only.sh      # AddressSanitizer
bash scripts/valgrind_core_only.sh  # Valgrind memcheck
```
