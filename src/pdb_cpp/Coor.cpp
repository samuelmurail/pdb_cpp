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

vector<string> Coor::get_uniq_chain_str() const {
    vector<string> chains;
    if (models_.empty()) {
        return chains;
    }
    vector<array<char, 2>> raw = models_[active_model].get_uniq_chain();
    chains.reserve(raw.size());
    for (const auto &chain : raw) {
        string chain_id;
        for (char letter : chain) {
            if (letter != '\0' && letter != ' ') {
                chain_id.push_back(letter);
            }
        }
        chains.push_back(chain_id);
    }
    return chains;
}

vector<string> Coor::get_aa_sequences(bool gap_in_seq, size_t frame) const {
    // Ensure the frame index is valid
    if (frame >= model_size()) {
        throw out_of_range("Frame index out of range");
    }
    if (models_[frame].size() == 0) {
        return {};
    }

    // Get the indexes of the selected atoms from the specified model
    vector<bool> CA_indexes = models_[frame].select_atoms("name CA");
    vector<array<char, 5>> resname_array = models_[frame].get_resname();
    vector<array<char, 2>> chain_array = models_[frame].get_chain();
    vector<int> resid_array = models_[frame].get_resid();

    if (resname_array.empty() || chain_array.empty() || resid_array.empty()) {
        return {};
    }

    auto chain_to_string = [](const array<char, 2> &chain) {
        string chain_id;
        for (char letter : chain) {
            if (letter != '\0' && letter != ' ') {
                chain_id.push_back(letter);
            }
        }
        return chain_id;
    };

    vector<string> seq_vec;
    unordered_map<string, size_t> chain_index;
    unordered_map<string, int> last_resid;
    size_t gap_num = 0;

    for (size_t i = 0; i < CA_indexes.size(); ++i) {
        if (!CA_indexes[i]) {
            continue;
        }

        string chain_id = chain_to_string(chain_array[i]);
        auto it = chain_index.find(chain_id);
        if (it == chain_index.end()) {
            chain_index[chain_id] = seq_vec.size();
            seq_vec.emplace_back("");
            last_resid[chain_id] = resid_array[i];
            it = chain_index.find(chain_id);
        }

        int prev_resid = last_resid[chain_id];
        if (gap_in_seq) {
            int diff = resid_array[i] - prev_resid;
            if (diff > 1) {
                gap_num = static_cast<size_t>(diff - 1);
                for (size_t j = 0; j < gap_num; ++j) {
                    seq_vec[it->second] += "-";
                }
            }
        }

        seq_vec[it->second] += convert_to_one_letter_resname(resname_array[i]);
        last_resid[chain_id] = resid_array[i];
    }

    return seq_vec;
}

vector<string> Coor::get_aa_sequences_dl(bool gap_in_seq, size_t frame) const {
    if (frame >= model_size()) {
        throw out_of_range("Frame index out of range");
    }
    if (models_[frame].size() == 0) {
        return {};
    }

    vector<bool> CA_indexes = models_[frame].select_atoms("name CA");
    vector<array<char, 5>> resname_array = models_[frame].get_resname();
    vector<array<char, 2>> chain_array = models_[frame].get_chain();
    vector<int> resid_array = models_[frame].get_resid();

    if (resname_array.empty() || chain_array.empty() || resid_array.empty()) {
        return {};
    }

    auto chain_to_string = [](const array<char, 2> &chain) {
        string chain_id;
        for (char letter : chain) {
            if (letter != '\0' && letter != ' ') {
                chain_id.push_back(letter);
            }
        }
        return chain_id;
    };

    vector<string> seq_vec;
    unordered_map<string, size_t> chain_index;
    unordered_map<string, int> last_resid;
    size_t gap_num = 0;

    for (size_t i = 0; i < CA_indexes.size(); ++i) {
        if (!CA_indexes[i]) {
            continue;
        }

        string chain_id = chain_to_string(chain_array[i]);
        auto it = chain_index.find(chain_id);
        if (it == chain_index.end()) {
            chain_index[chain_id] = seq_vec.size();
            seq_vec.emplace_back("");
            last_resid[chain_id] = resid_array[i];
            it = chain_index.find(chain_id);
        }

        int prev_resid = last_resid[chain_id];
        if (gap_in_seq) {
            int diff = resid_array[i] - prev_resid;
            if (diff > 1) {
                gap_num = static_cast<size_t>(diff - 1);
                for (size_t j = 0; j < gap_num; ++j) {
                    seq_vec[it->second] += "-";
                }
            }
        }

        seq_vec[it->second] += convert_to_one_letter_resname_dl(resname_array[i]);
        last_resid[chain_id] = resid_array[i];
    }

    return seq_vec;
}