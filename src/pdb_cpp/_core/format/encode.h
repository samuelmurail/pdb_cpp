#ifndef PDB_CPP_FORMAT_ENCODE_H
#define PDB_CPP_FORMAT_ENCODE_H

/**
 * @file encode.h
 * @brief Hybrid-36 encoding helpers for PDB-style fields.
 */

#include <string>

/** Decode a hybrid-36 encoded field. */
int hy36decode(int width, const std::string& s);
/** Encode an integer using hybrid-36. */
std::string hy36encode(int width, int value);

#endif // PDB_CPP_FORMAT_ENCODE_H
