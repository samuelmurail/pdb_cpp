#ifndef PDB_CPP_FORMAT_ENCODE_H
#define PDB_CPP_FORMAT_ENCODE_H

#include <string>

int hy36decode(int width, const std::string& s);
std::string hy36encode(int width, int value);

#endif // PDB_CPP_FORMAT_ENCODE_H
