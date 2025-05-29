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
print('pdb_numpy alignement')
alignement.print_align_seq(align[0], align[1], line_len=80)

print('pdb_cpp alignement')

align_seq_1, align_seq_2, _ = alignement.align_seq(seq_1, seq_2)
alignement.print_align_seq(align_seq_1, align_seq_2, line_len=80)


seq_1 = (
    "AQDMVSPPPPIADEPLTVNTGIYLIECYSLDDKAETFKVNAFLSLSWKDRRLAFDPV"
    "RSGVRVKTYEPEAIWIPEIRFVNVENARDADVVDISVSPDGTVQYLERFSARVLSPLDFRRYPFDSQTLHIYLIVRSV"
    "DTRNIVLAVDLEKVGKNDDVFLTGWDIESFTAVVKPANFALEDRLESKLDYQLRISRQYFSYIPNIILPMLFILFISW"
    "TAFWSTSYEANVTLVVSTLIAHIAFNILVETNLPKTPYMTYTGAIIFMIYLFYFVAVIEVTVQHYLKVESQPARAASI"
    "TRASRIAFPVVFLLANIILAFLFFGF"
)
seq_2 = (
    "APSEFLDKLMGKVSGYDARIRPNFKGPPVNVTCNIFINSFGSIAETTMDYRVNIFLR"
    "QQWNDPRLAYSEYPDDSLDLDPSMLDSIWKPDLFFANEKGANFHEVTTDNKLLRISKNGNVLYSIRITLVLACPMDLK"
    "NFPMDVQTCIMQLESFGYTMNDLIFEWDEKGAVQVADGLTLPQFILKEEKDLRYCTKHYNTGKFTCIEARFHLERQMG"
    "YYLIQMYIPSLLIVILSWVSFWINMDAAPARVGLGITTVLTMTTQSSGSRASLPKVSYVKAIDIWMAVCLLFVFSALL"
    "EYAAVNFIARAGTKLFISRAKRIDTVSRVAFPLVFLIFNIFYWITYKLVPR"
)

seq_3 = "VSPPPPIADEPLTVNTGIYLIECYSLDDKAETFKVNAFLSLSWKDRRLAFDPVRSGVRVKTYEPEAIWIPEIRFVNVENARDADVVDISVSPDGTVQYLERFSARVLSPLDFRRYPFDSQTLHIYLIVRSVDTRNIVLAVDLEKVGKNDDVFLTGWDIESFTAVVKPANFALEDRLESKLDYQLRISRQYFSYIPNIILPMLFILFISWTAFWSTSYEANVTLVVSTLIAHIAFNILVETNLPKTPYMTYTGAIIFMIYLFYFVAVIEVTVQHYLKVESQPARAASITRASRIAFPVVFLLANIILAFLFF"
seq_4 = "LSPSDFLDKLMGRTSGYDARIRPNFKGPPVNVTCNIFINSFGSVTETTMDYRVNIFLRQQWNDSRLAYSEYPDDSLDLDPSMLDSIWKPDLFFANEKGANFHDVTTDNKLLRISKNGKVLYSIRLTLTLSCPMDLKNFPMDVQTCTMQLESFGYTMNDLIFEWLSDGPVQVAEGLTLPQFILKEEKELGYCTKHYNTGKFTCIEVKFHLERQMGYYLIQMYIPSLLIVILSWVSFWINMDAAPARVALGITTVLTMTTQSSGSRASLPKVSYVKAIDIWMAVCLLFVFAALLEYAAVNFVSRKFVDRAKRIDTISRAAFPLAFLIFNIFYWITYKIIRG"

import pdb_numpy
import pdb_numpy.alignement
align = pdb_numpy.alignement.align_seq_cython(seq_3, seq_4)
print('pdb_numpy alignement')
alignement.print_align_seq(align[0], align[1], line_len=80)

print('pdb_cpp alignement')

align_seq_1, align_seq_2, score = alignement.align_seq(seq_3, seq_4)
print(f"Score: {score}")
alignement.print_align_seq(align_seq_1, align_seq_2, line_len=80)

file_name = "3eam.pdb"
coor1 = Coor(file_name)

file_name = "5bkg.pdb"
coor2 = Coor(file_name)

seq_1 = coor1.get_aa_seq()
seq_2 = coor2.get_aa_seq()

align_seq_3, align_seq_4, score = alignement.align_seq(seq_1['A'], seq_2['A'])
print(f"Score: {score}")
alignement.print_align_seq(align_seq_3, align_seq_4, line_len=80)

print()

test = alignement.get_common_atoms(
    coor_1,
    coor_2);
