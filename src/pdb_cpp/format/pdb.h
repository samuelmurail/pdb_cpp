#ifndef PDB_CPP_FORMAT_PDB_H
#define PDB_CPP_FORMAT_PDB_H

#include <string>
#include "../Coor.h"

class Coor;

Coor PDB_parse(const string& filename);
bool PDB_write(const Coor& coor, const string& filename);

#endif // PDB_CPP_FORMAT_PDB_H
