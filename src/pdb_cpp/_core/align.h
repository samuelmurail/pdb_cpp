#ifndef SEQ_ALIGN_H
#define SEQ_ALIGN_H

/**
 * @file align.h
 * @brief Sequence-alignment and coordinate-alignment helpers.
 */

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

/** Align two coordinate structures using precomputed atom index lists. */
void coor_align(Coor& coor_1, const Coor& coor_2, 
               const std::vector<int>& index_1, 
               const std::vector<int>& index_2, 
               int frame_ref = 0);

/** Compute RMSD for two atom-index selections. */
float rmsd(
    const Model &model_1,
    const Model &model_2,
    const std::vector<int>& index_1, 
    const std::vector<int>& index_2);

/** Align chains by sequence and return the best RMSD plus index mappings. */
std::tuple<std::vector<float>, std::vector<int>, std::vector<int>> align_seq_based(
    Coor &coor_1,
    const Coor &coor_2,
    const std::vector<std::string> &chain_1= {"A"},
    const std::vector<std::string> &chain_2= {"A"},
    const std::vector<std::string> &back_names= {"C", "N", "O", "CA"},
    const std::string &matrix_file="",
    const int frame_ref=0
);

/** Align structures using explicit atom index lists. */
std::tuple<std::vector<float>, std::vector<int>, std::vector<int>> align_index_based(
    Coor &coor_1,
    const Coor &coor_2,
    const std::vector<int> &index_1,
    const std::vector<int> &index_2,
    const int frame_ref=0
);

/** Try all chain permutations and return the best-scoring mapping. */
std::pair<std::vector<float>, std::pair<std::vector<int>, std::vector<int>>> align_chain_permutation(
    const Coor &coor_1,
    const Coor &coor_2,
    const std::vector<std::string> &back_names= {"C", "N", "O", "CA"},
    const std::string &matrix_file="",
    const int frame_ref=0
);



#endif // SEQ_ALIGN_H