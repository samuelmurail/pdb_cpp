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
    active_model = 0;
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
    selected.active_model = active_model;

    for (size_t i = 0; i < models_.size(); ++i) {
        Model model = models_[i].select_index(indexes);
        selected.add_Model(model);
    }

    return selected;
}

vector<array<char, 2>> Coor::get_uniq_chain() const {
    if (models_.empty()) {
        throw runtime_error("No models available");
    }
    return models_[active_model].get_uniq_chain();
}

vector<string> Coor::get_aa_sequences(bool gap_in_seq, size_t frame) const {
    // Ensure the frame index is valid
    if (frame >= model_size()) {
        throw out_of_range("Frame index out of range");
    }

    // Get the indexes of the selected atoms from the specified model
    vector<bool> CA_indexes = models_[frame].select_atoms("name CA");
    vector<array<char, 5>> resname_array = models_[frame].get_resname();
    vector<array<char, 2>> chain_array = models_[frame].get_chain();
    vector<int> resid_array = models_[frame].get_resid();

    vector<array<char, 2>> uniq_chains= get_uniq_chain();

    array<char, 2> old_chain = chain_array[0];
    // Get the index of the old chain in the unique chains
    auto it = find(uniq_chains.begin(), uniq_chains.end(), old_chain);
    if (it == uniq_chains.end()) {
        throw runtime_error("Chain not found in unique chains");
    }
    int chain_index = distance(uniq_chains.begin(), it);

    vector<string> seq_vec;
    int old_resid = resid_array[0];
    seq_vec.emplace_back("");

    for (size_t i = 0; i < CA_indexes.size(); ++i) {
        if (CA_indexes[i]) {
            if (chain_array[i] != old_chain) {
                // New chain or new residue
                old_chain = chain_array[i];
                old_resid = resid_array[i];
                // Get the index of the old chain in the unique chains
                it = find(uniq_chains.begin(), uniq_chains.end(), old_chain);
                if (it == uniq_chains.end()) {
                    throw runtime_error("Chain not found in unique chains");
                }
                chain_index = distance(uniq_chains.begin(), it);
                seq_vec.emplace_back("");
            }
            if (resid_array[i] != old_resid) {
                // New residue
                old_resid = resid_array[i];
                if (gap_in_seq) {
                    seq_vec[chain_index] += "-"; // Add gap for new residue
                }
            }
            seq_vec[chain_index] += convert_to_one_letter_resname(resname_array[i]);
            old_resid += 1;
        }
    }
    return seq_vec;
}