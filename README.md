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
