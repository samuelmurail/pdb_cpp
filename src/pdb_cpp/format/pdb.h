#ifndef PDB_CPP_FORMAT_PDB_H
#define PDB_CPP_FORMAT_PDB_H

#include <string>
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include "../Coor.h"
#include "../Model.h"

class Coor;

Coor PDB_parse(const std::string& filename);
bool PDB_write(const Coor& coor, const std::string& filename);

#endif // PDB_CPP_FORMAT_PDB_H
