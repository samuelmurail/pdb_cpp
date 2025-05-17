#ifndef PDB_CPP_FORMAT_PDB_H
#define PDB_CPP_FORMAT_PDB_H

#include <string>
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include "../Model.h"

class Model;

Model PDB_parse(const std::string& filename);

#endif // PDB_CPP_FORMAT_PDB_H
