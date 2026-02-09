#ifndef SEQ_ALIGN_H
#define SEQ_ALIGN_H

#include <vector>
#include <string>
#include <tuple>
#include <utility>
#include <cassert>

#include "Coor.h"

std::pair<std::vector<int>, std::vector<int>> get_common_atoms(
    const Coor &coor_1,
    const Coor &coor_2,
    const std::vector<std::string> &chain_1= {"A"},
    const std::vector<std::string> &chain_2= {"A"},
    const std::vector<std::string> &back_names= {"C", "N", "O", "CA"},
    const std::string &matrix_file=""
);

// Function to align two coordinate structures
void coor_align(Coor& coor_1, const Coor& coor_2, 
               const std::vector<int>& index_1, 
               const std::vector<int>& index_2, 
               int frame_ref = 0);

float rmsd(
    const Model &model_1,
    const Model &model_2,
    const std::vector<int>& index_1, 
    const std::vector<int>& index_2);

std::tuple<std::vector<float>, std::vector<int>, std::vector<int>> align_seq_based(
    Coor &coor_1,
    const Coor &coor_2,
    const std::vector<std::string> &chain_1= {"A"},
    const std::vector<std::string> &chain_2= {"A"},
    const std::vector<std::string> &back_names= {"C", "N", "O", "CA"},
    const std::string &matrix_file="",
    const int frame_ref=0
);

std::pair<std::vector<float>, std::pair<std::vector<int>, std::vector<int>>> align_chain_permutation(
    const Coor &coor_1,
    const Coor &coor_2,
    const std::vector<std::string> &back_names= {"C", "N", "O", "CA"},
    const std::string &matrix_file="src/pdb_cpp/data/blosum62.txt",
    const int frame_ref=0
);



#endif // SEQ_ALIGN_H