#include <cstring>
#include <iomanip>
#include <unordered_map>

#include "Coor.h"
#include "Model.h"
#include "format/pdb.h"
#include "sequence.h"

using namespace std;

bool endswith(string const &fullString, string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

void Coor::clear()
{
    models_.clear();
    crystal_pack.clear();
    active_model_ = 0;
}

bool Coor::read(const string &filename) {

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

bool Coor::write(const string &filename) const {
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

Coor Coor::select_atoms(const string &selection, size_t frame) const {
    // Ensure the frame index is valid
    if (frame >= model_size()) {
        throw out_of_range("Frame index out of range");
    }

    // Get the indexes of the selected atoms from the specified model
    vector<bool> indexes = models_[frame].select_atoms(selection);

    // Return a new Coor object with the selected atoms
    return select_bool_index(indexes);
}

Coor Coor::select_bool_index(const vector<bool> &indexes) const {
    Coor selected;
    selected.clear();
    selected.crystal_pack = crystal_pack;
    selected.transformation = transformation;
    selected.symmetry = symmetry;
    selected.active_model_ = active_model_;

    for (size_t i = 0; i < models_.size(); ++i) {
        Model model = models_[i].select_index(indexes);
        selected.add_Model(model);
    }

    return selected;
}

void Coor::get_aa_seq(bool gap_in_seq, size_t frame) const {
    // Ensure the frame index is valid
    if (frame >= model_size()) {
        throw out_of_range("Frame index out of range");
    }

    // Get the indexes of the selected atoms from the specified model
    vector<bool> CA_indexes = models_[frame].select_atoms("name CA");
    vector<array<char, 5>> resname_array = models_[frame].get_resname();

    unordered_map<char, char> AA_DICT_L = {


    for (size_t i = 0; i < CA_indexes.size(); ++i) {
        if (CA_indexes[i]) {
            cout << convert_to_one_letter_resname(resname_array[i]);
        }
    }
}