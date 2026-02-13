# pdb_cpp

Library to use pdb/mmcif files with c++.


## Installation

```bash
git clone https://github.com/samuelmurail/pdb_cpp
cd pdb_cpp
python -m pip install -e .
```

For development checks:

```bash
python -m pip install -r requirements.txt
pytest
```

## Usage
```python
import pdb_cpp

pdb_cpp.read_pdb("1aon.pdb")
```


## Adding c++ code

To add c++ code, create new files in `src/pdb_cpp/_core` and then modify `setup.py` to include them in the build process.

```
ext_modules = [
    Pybind11Extension(
        "core",
        [
            "src/pdb_cpp/_core/pybind.cpp",
            "src/pdb_cpp/_core/sequence_align.cpp",  # Include sequence_align.cpp
            "src/pdb_cpp/_core/Model.cpp",
            "src/pdb_cpp/_core/Coor.cpp",
            # Add other source files as needed
        ],
        include_dirs=["src/pdb_cpp/_core"],  # Include directory for headers
        cxx_std=17,  # Use C++17 standard
    ),
]
```

Also add the function and class declarations in the `src/pdb_cpp/_core/pybind.cpp` file to expose your new functionality to Python.
