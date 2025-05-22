import time
import pdb_numpy
import sys

sys.path.insert(0, "./src")

import pdb_cpp

test = pdb_cpp.Coor("3eam.pdb")
print(test.len)