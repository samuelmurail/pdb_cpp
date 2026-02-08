#include <cstring>
#include <iomanip>
#include <unordered_map>
#include <string>
#include <fstream>

#include "Coor.h"
#include "Model.h"
#include "format/pdb.h"
#include "format/mmcif.h"
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
    conect.clear();
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
    if (endswith(filename, ".cif")) {
        *this = MMCIF_parse(filename);
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
    if (endswith(filename, ".cif")) {
        MMCIF_write(*this, filename);
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

    unordered_map<int, int> index_dict;
    const vector<int> &num_list = models_[active_model].get_num();
    int new_serial = 1;
    for (size_t i = 0; i < indexes.size(); ++i) {
        if (indexes[i]) {
            index_dict[num_list[i]] = new_serial;
            ++new_serial;
        }
    }

    for (size_t i = 0; i < models_.size(); ++i) {
        Model model = models_[i].select_index(indexes);
        selected.add_Model(model);
    }

    for (size_t m = 0; m < selected.models_.size(); ++m) {
        for (size_t i = 0; i < selected.models_[m].size(); ++i) {
            selected.models_[m].set_num(i, static_cast<int>(i + 1));
        }
    }

    selected.conect.clear();
    for (const auto &kv : conect) {
        auto key_it = index_dict.find(kv.first);
        if (key_it == index_dict.end()) {
            continue;
        }
        vector<int> new_values;
        new_values.reserve(kv.second.size());
        for (int val : kv.second) {
            auto val_it = index_dict.find(val);
            if (val_it != index_dict.end()) {
                new_values.push_back(val_it->second);
            }
        }
        if (!new_values.empty()) {
            selected.conect[key_it->second] = move(new_values);
        }
    }

    return selected;
}

vector<int> Coor::get_index_select(const string selection, size_t frame) const{
    // Ensure the frame index is valid
    if (frame >= model_size()) {
        throw out_of_range("Frame index out of range");
    }

    // Get the indexes of the selected atoms from the specified model
    return models_[frame].get_index_select(selection);
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
    size_t gap_num;

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
                if (gap_in_seq) {
                    gap_num = resid_array[i] - old_resid;
                    for (size_t j = 0; j < gap_num; ++j) {
                        seq_vec[chain_index] += "-"; // Add gap for new residue
                    }
                }
                old_resid = resid_array[i];
            }
            seq_vec[chain_index] += convert_to_one_letter_resname(resname_array[i]);
            old_resid += 1;
        }
    }
    return seq_vec;
}