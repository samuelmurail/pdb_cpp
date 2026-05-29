#ifndef PDB_CPP_FORMAT_PDB_H
#define PDB_CPP_FORMAT_PDB_H

/**
 * @file pdb.h
 * @brief PDB reader and writer entry points.
 */

#include <string>
#include "../Coor.h"

class Coor;

/** Parse a PDB file into a Coor object. */
Coor PDB_parse(const std::string& filename);
/** Serialize a Coor object to a PDB file. */
bool PDB_write(const Coor& coor, const std::string& filename);

#endif // PDB_CPP_FORMAT_PDB_H
