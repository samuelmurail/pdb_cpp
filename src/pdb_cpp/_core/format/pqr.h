#ifndef PDB_CPP_FORMAT_PQR_H
#define PDB_CPP_FORMAT_PQR_H

#include <string>
#include "../Coor.h"

class Coor;

Coor PQR_parse(const std::string& filename);
bool PQR_write(const Coor& coor, const std::string& filename);

#endif // PDB_CPP_FORMAT_PQR_H
