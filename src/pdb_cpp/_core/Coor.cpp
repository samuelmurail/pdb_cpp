#include <cstring>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <fstream>
#include <cctype>

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
        return model_size() > 0;
    }
    if (endswith(filename, ".cif")) {
        *this = MMCIF_parse(filename);
        return model_size() > 0;
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
            selected.conect[key_it->second] = std::move(new_values);
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

    vector<int> CA_indexes = models_[frame].get_index_select("name CA and not altloc B C D E F");
    vector<array<char, 5>> resname_array = models_[frame].get_resname();
    vector<array<char, 2>> chain_array = models_[frame].get_chain();
    vector<int> resid_array = models_[frame].get_resid();
    vector<int> uniq_resid_array = models_[frame].get_uniqresid();

    Coor n_c_cb = select_atoms("name N C CB and not altloc B C D E F", frame);

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

    for (int idx : CA_indexes) {
        string chain_id = chain_to_string(chain_array[idx]);
        auto it = chain_index.find(chain_id);
        if (it == chain_index.end()) {
            chain_index[chain_id] = seq_vec.size();
            seq_vec.emplace_back("");
            last_resid[chain_id] = resid_array[idx];
            it = chain_index.find(chain_id);
        }

        int prev_resid = last_resid[chain_id];
        if (gap_in_seq) {
            int diff = resid_array[idx] - prev_resid;
            if (diff > 1) {
                gap_num = static_cast<size_t>(diff - 1);
                for (size_t j = 0; j < gap_num; ++j) {
                    seq_vec[it->second] += "-";
                }
            }
        }

        char aa = convert_to_one_letter_resname_any(resname_array[idx]);
        if (aa == 'G') {
            seq_vec[it->second] += aa;
        } else {
            int uniq_resid = uniq_resid_array[idx];
            vector<int> n_index = n_c_cb.get_index_select(
                "name N and residue " + to_string(uniq_resid), frame);
            vector<int> c_index = n_c_cb.get_index_select(
                "name C and residue " + to_string(uniq_resid), frame);
            vector<int> cb_index = n_c_cb.get_index_select(
                "name CB and residue " + to_string(uniq_resid), frame);

            if (!n_index.empty() && !c_index.empty() && !cb_index.empty()) {
                const Model ncc_model = n_c_cb.get_Models(static_cast<int>(frame));
                auto to_array = [](const Model &model, int index) {
                    return std::array<float, 3>{
                        model.get_x()[index],
                        model.get_y()[index],
                        model.get_z()[index]
                    };
                };

                const Model ref_model = get_Models(static_cast<int>(frame));
                std::array<float, 3> ca_xyz = {
                    ref_model.get_x()[idx],
                    ref_model.get_y()[idx],
                    ref_model.get_z()[idx]
                };

                float dihed = atom_dihed_angle(
                    ca_xyz,
                    to_array(ncc_model, n_index[0]),
                    to_array(ncc_model, c_index[0]),
                    to_array(ncc_model, cb_index[0]));

                if (dihed > 0.0f) {
                    seq_vec[it->second] += aa;
                } else {
                    seq_vec[it->second] += static_cast<char>(
                        std::tolower(static_cast<unsigned char>(aa)));
                }
            } else {
                seq_vec[it->second] += aa;
            }
        }
        last_resid[chain_id] = resid_array[idx];
    }

    return seq_vec;
}

unordered_map<string, string> Coor::get_aa_seq(bool gap_in_seq, size_t frame) const {
    vector<string> seq_vec = get_aa_sequences(gap_in_seq, frame);
    if (seq_vec.empty()) {
        return {};
    }

    vector<int> ca_index = models_[frame].get_index_select("name CA");
    vector<array<char, 2>> chain_array = models_[frame].get_chain();

    auto chain_to_string = [](const array<char, 2> &chain) {
        string chain_id;
        for (char letter : chain) {
            if (letter != '\0' && letter != ' ') {
                chain_id.push_back(letter);
            }
        }
        return chain_id;
    };

    vector<string> chain_keys;
    unordered_set<string> seen;
    for (int idx : ca_index) {
        string chain_id = chain_to_string(chain_array[idx]);
        if (seen.insert(chain_id).second) {
            chain_keys.push_back(chain_id);
        }
    }

    unordered_map<string, string> seq_dict;
    size_t count = min(chain_keys.size(), seq_vec.size());
    for (size_t i = 0; i < count; ++i) {
        seq_dict[chain_keys[i]] = seq_vec[i];
    }
    return seq_dict;
}

unordered_map<string, string> Coor::get_aa_DL_seq(bool gap_in_seq, size_t frame) const {
    vector<string> seq_vec = get_aa_sequences_dl(gap_in_seq, frame);
    if (seq_vec.empty()) {
        return {};
    }

    vector<int> ca_index = models_[frame].get_index_select("name CA and not altloc B C D E F");
    vector<array<char, 2>> chain_array = models_[frame].get_chain();

    auto chain_to_string = [](const array<char, 2> &chain) {
        string chain_id;
        for (char letter : chain) {
            if (letter != '\0' && letter != ' ') {
                chain_id.push_back(letter);
            }
        }
        return chain_id;
    };

    vector<string> chain_keys;
    unordered_set<string> seen;
    for (int idx : ca_index) {
        string chain_id = chain_to_string(chain_array[idx]);
        if (seen.insert(chain_id).second) {
            chain_keys.push_back(chain_id);
        }
    }

    unordered_map<string, string> seq_dict;
    size_t count = min(chain_keys.size(), seq_vec.size());
    for (size_t i = 0; i < count; ++i) {
        seq_dict[chain_keys[i]] = seq_vec[i];
    }
    return seq_dict;
}

unordered_map<string, string> Coor::get_aa_na_seq(bool gap_in_seq, size_t frame) const {
    if (frame >= model_size()) {
        throw out_of_range("Frame index out of range");
    }
    if (models_[frame].size() == 0) {
        return {};
    }

    vector<bool> sel = models_[frame].select_atoms(
        "((protein and name CA) or (dna and name P)) and not altloc B C D E F");
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

    unordered_map<string, string> seq_dict;
    unordered_map<string, int> last_resid;

    for (size_t i = 0; i < sel.size(); ++i) {
        if (!sel[i]) {
            continue;
        }

        string chain_id = chain_to_string(chain_array[i]);
        if (seq_dict.find(chain_id) == seq_dict.end()) {
            seq_dict[chain_id] = "";
            last_resid[chain_id] = resid_array[i];
        }

        if (gap_in_seq) {
            int diff = resid_array[i] - last_resid[chain_id];
            if (diff > 1 && !seq_dict[chain_id].empty()) {
                seq_dict[chain_id].append(static_cast<size_t>(diff - 1), '-');
            }
        }

        try {
            seq_dict[chain_id].push_back(convert_to_one_letter_resname_na(resname_array[i]));
        } catch (const std::exception &) {
            seq_dict[chain_id].push_back('X');
        }

        last_resid[chain_id] = resid_array[i];
    }

    return seq_dict;
}

Coor Coor::remove_incomplete_backbone_residues(const vector<string> &back_atom) const {
    Coor no_alter_loc = select_atoms("protein and not altloc B C D");
    if (no_alter_loc.model_size() == 0) {
        return no_alter_loc;
    }

    string name_sel = "name";
    for (const auto &atom : back_atom) {
        name_sel += " " + atom;
    }
    Coor backbone = no_alter_loc.select_atoms(name_sel);
    if (backbone.model_size() == 0) {
        return no_alter_loc;
    }

    Model bb_model = backbone.get_Models(0);
    const vector<int> &uniq_resid = bb_model.get_uniqresid();
    if (uniq_resid.empty()) {
        return no_alter_loc;
    }

    unordered_map<int, int> counts;
    for (int resid : uniq_resid) {
        counts[resid] += 1;
    }

    vector<int> to_remove;
    to_remove.reserve(counts.size());
    int expected = static_cast<int>(back_atom.size());
    for (const auto &kv : counts) {
        if (kv.second != expected) {
            to_remove.push_back(kv.first);
        }
    }

    if (to_remove.empty()) {
        return no_alter_loc;
    }

    string remove_sel = "not residue";
    for (int resid : to_remove) {
        remove_sel += " " + to_string(resid);
    }
    return no_alter_loc.select_atoms(remove_sel);
}