import time
import pdb_numpy
import sys

sys.path.insert(0, "./src")

from pdb_cpp import Coor

file_name = "3eam.pdb"
coor = Coor(file_name)
seqs = coor.get_aa_seq()
print(seqs)