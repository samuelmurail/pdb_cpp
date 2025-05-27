#ifndef SEQ_ALIGN_H
#define SEQ_ALIGN_H

#include <vector>
#include <string>

struct Alignment {
    std::string seq1;
    std::string seq2;
    int score;
};


Alignment sequence_align(
    const std::string seq1,
    const std::string seq2,
    const std::string matrix_file,
    int GAP_COST=-11,
    int GAP_EXT=+1);

#endif // SEQ_ALIGN_H