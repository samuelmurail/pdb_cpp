#ifndef SEQ_ALIGN_H
#define SEQ_ALIGN_H

#include <vector>
#include <string>
#include <utility>
#include <cassert>


#include "Coor.h"

struct Alignment {
    std::string seq1;
    std::string seq2;
    int score;
};

void print_alignment(Alignment &alignment);

Alignment sequence_align(
    const std::string seq1,
    const std::string seq2,
    const std::string matrix_file ="",
    int GAP_COST=-11,
    int GAP_EXT=-1);

// pair<vector<int>, vector<int>> get_common_atoms(
//     const Coor &coor_1,
//     const Coor &coor_2,
//     const std::vector<std::string> &chain_1 = {"A"},
//     const std::vector<std::string> &chain_2 = {"A"},
//     const std::vector<std::string> &back_names = {"C", "N", "O", "CA"},
//     const std::string matrix_file ="");


//pair<vector<int>, vector<int>> get_common_atoms(
std::pair<std::vector<int>, std::vector<int>> get_common_atoms(
    const Coor &coor_1,
    const Coor &coor_2,
    const std::vector<std::string> &chain_1= {"A"},
    const std::vector<std::string> &chain_2= {"A"},
    const std::vector<std::string> &back_names= {"C", "N", "O", "CA"},
    const std::string &matrix_file=""
);

// Function to align two coordinate structures
void coor_align(Coor& coor_1, Coor& coor_2, 
               const std::vector<int>& index_1, 
               const std::vector<int>& index_2, 
               int frame_ref = 0);

    
#endif // SEQ_ALIGN_H