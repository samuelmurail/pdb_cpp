#include <cstring>
#include <iomanip>

#include "Coor.h"
#include "format/pdb.h"

using namespace std;


bool endswith (string const &fullString, string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}


void Coor::clear() {
    models_.clear();
    crystal_pack.clear();
    active_model_ = 0;
}


bool Coor::read(const string& filename) {
    
    ifstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << endl;
        return false;
    }
    clear();

    if (endswith(filename, ".pdb")) {
        // cout << "Reading PDB file: " << filename << endl;
        *this = PDB_parse(filename);
        return true;
    } 
    return false;

}

bool Coor::write(const string& filename) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << endl;
        return false;
    }

    if (endswith(filename, ".pdb")) {
        // cout << "Writing PDB file: " << filename << endl;
        PDB_write(*this, filename);
        return true;
    } 
    return false;
    
}