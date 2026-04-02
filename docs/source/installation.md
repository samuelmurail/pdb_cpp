# Installation

## Requirements

- **Python** >= 3.10
- **C++ compiler** with C++17 support (GCC >= 7, Clang >= 5, MSVC >= 2017)
- **NumPy** >= 1.22

Build-time only (installed automatically):

- **pybind11** == 2.13

## From PyPI (recommended)

```bash
pip install pdb-cpp
```

This installs a pre-built wheel when available for your platform, or
compiles from source otherwise.

## From source

```bash
git clone https://github.com/samuelmurail/pdb_cpp.git
cd pdb_cpp
pip install .
```

### Development install

For active development with editable installs:

```bash
pip install -r requirements.txt   # build dependencies
pip install -e .                   # editable install
```

After changing C++ code, rebuild the extension:

```bash
pip install -e . --no-build-isolation
```

## Verify the installation

```bash
python -c "from pdb_cpp import Coor; print('pdb_cpp imported successfully')"
```

## Run the test suite

```bash
pip install pytest
pytest
```

## Platform notes

### Linux

A C++ compiler is usually pre-installed. If not:

```bash
sudo apt-get install g++    # Debian/Ubuntu
sudo dnf install gcc-c++    # Fedora
```

### macOS

The Xcode command-line tools provide `clang++`:

```bash
xcode-select --install
```

### Windows

Install the [Visual Studio Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
and select the "C++ build tools" workload.