#ifndef PDB_CPP_FORMAT_GRO_H
#define PDB_CPP_FORMAT_GRO_H

#include <string>
#include "../Coor.h"

class Coor;

Coor GRO_parse(const std::string& filename);
bool GRO_write(const Coor& coor, const std::string& filename);

#endif // PDB_CPP_FORMAT_GRO_H
