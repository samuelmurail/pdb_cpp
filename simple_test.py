import time
import pdb_numpy
import sys

sys.path.insert(0, "./src")

from pdb_cpp import Coor, TMalign

file_name = "3eam.pdb"
coor = Coor(file_name)
seqs = coor.get_aa_seq()
print(seqs)

test = TMalign.compute_secondary_structure(coor)[0]

for key in test:
    print(f"{key}: {len(test[key])} {len(seqs[key])}")

print(test)