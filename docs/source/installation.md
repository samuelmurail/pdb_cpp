# Installation Quick Start

## Through Pypi

Install from PyPI:

```
pip install pdb-cpp
```

## Get sources from the GithubRepo

The sources for pdb_cpp can be downloaded from GitHub.

You can either clone the public repository:

```bash
$ git clone git@github.com:samuelmurail/pdb_cpp.git
```

Or download the tarball:

```bash
$ curl -OJL https://github.com/samuelmurail/pdb_cpp/tarball/master
```

Once you have a copy of the source, switch to the `pdb_cpp` directory.

```bash
$ cd pdb_cpp
```

##  Install `pdb_cpp`

Once you have a copy of the source and have created an environment, install with:

```bash
$ pip install .
```

For development mode:

```bash
$ pip install -e .
```

## Test Installation

To test the installation, use `pytest`:

```bash
$ pytest
```

If you changed C++ extension code, reinstall to rebuild the extension:

```bash
$ pip install -e . --no-build-isolation
```

and then run the tests again.