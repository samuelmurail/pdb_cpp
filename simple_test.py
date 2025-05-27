import time
import pdb_numpy
import sys

sys.path.insert(0, "./src")

from pdb_cpp import Coor, TMalign, alignement

file_name = "3eam.pdb"
coor = Coor(file_name)
seqs = coor.get_aa_seq()
print(seqs)

test = TMalign.compute_secondary_structure(coor)[0]

for key in test:
    print(f"{key}: {len(test[key])} {len(seqs[key])}")

print(test)

seq_1 = "AQDMVSPPXPIADEPLTVXSLSWKDRRL"
seq_2 = "AQDMVSPPPPIADEPLTVNTGIYLIECYSLDDKAETFKVNAFLSLSWKDRRLAFDPV"

import pdb_numpy
import pdb_numpy.alignement
align = pdb_numpy.alignement.align_seq_cython(seq_1, seq_2)
print(align)

align_seq_1, align_seq_2 = alignement.align_seq(seq_1, seq_2)
print(align_seq_1, align_seq_2)
alignement.print_align_seq(align_seq_1, align_seq_2, line_len=80)
