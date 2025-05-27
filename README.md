# pdb_cpp

Library to use pdb/mmcif files with c++.


## Installation

```bash
git clone https://github.com/samuelmurail/pdb_cpp
cd pdb_cpp
python setup.py build_ext --inplace
```

## Usage
```python
import pdb_cpp

pdb_cpp.read_pdb("1aon.pdb")
```


## Adding c++ code

To add c++ code, you can create a new file in the `src` directory and add your code there. Then, you can modify the `setup.py` file to include your new file in the build process.

```
ext_modules = [
    Pybind11Extension(
        "core",
        [
            "src/pdb_cpp/pybind.cpp",
            "src/pdb_cpp/sequence_align.cpp",  # Include sequence_align.cpp
            "src/pdb_cpp/Model.cpp",
            "src/pdb_cpp/Coor.cpp",
            # Add other source files as needed
        ],
        include_dirs=["src/pdb_cpp"],  # Include directory for headers
        cxx_std=17,  # Use C++17 standard
    ),
]
```

Also add the function and class declarations in the `src/pdb_cpp/pybind.cpp` file to expose your new functionality to Python.
