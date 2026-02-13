#ifndef PDB_CPP_FORMAT_MMCIF_H
#define PDB_CPP_FORMAT_MMCIF_H

#include <string>
#include "../Coor.h"

class Coor;

Coor MMCIF_parse(const std::string &filename);
std::string get_mmcif_string(const Coor &coor);
bool MMCIF_write(const Coor &coor, const std::string &filename);

#endif // PDB_CPP_FORMAT_MMCIF_H
