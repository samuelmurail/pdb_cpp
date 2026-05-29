#ifndef SEQUENCE_H
#define SEQUENCE_H

/**
 * @file sequence.h
 * @brief Residue-to-sequence conversion helpers.
 */

#include <array>
#include <cstring>
#include <stdexcept>

/** Convert a residue name to a one-letter amino-acid code. */
char convert_to_one_letter_resname(const std::array<char, 5> &resname_array);
/** Convert a residue name to a one-letter code, preserving D-residue case. */
char convert_to_one_letter_resname_dl(const std::array<char, 5> &resname_array);
/** Convert a residue name to a one-letter nucleic-acid code. */
char convert_to_one_letter_resname_na(const std::array<char, 5> &resname_array);
/** Convert a residue name to the most appropriate one-letter code. */
char convert_to_one_letter_resname_any(const std::array<char, 5> &resname_array);

#endif // SEQUENCE_H