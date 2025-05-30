#include <cstring>
#include <iomanip>

#include "Model.h"
#include "select.h"
#include "geom.h"

using namespace std;

void Model::clear() {
    x_.clear();
    y_.clear();
    z_.clear();
    name_.clear();
    resname_.clear();
    resid_.clear();
    chain_.clear();
    occ_.clear();
    beta_.clear();
    alterloc_.clear();
    insertres_.clear();
    elem_.clear();
    num_.clear();
    field_.clear();
    uniqresid_.clear();
}

size_t Model::size() const {
    return x_.size();
}

bool Model::addAtom(
    int num,
    const array<char, 5> &name_array,
    const array<char, 5> &resname_array,
    int res_id,
    const array<char, 2> &chain_array,
    float x, float y, float z, float occ, float beta,
    const array<char, 2> &alterloc = {' ', '\0'},
    const array<char, 5> &elem = {' ', ' ', ' ', ' ', '\0'},
    const array<char, 2> &insertres = {' ', '\0'},
    bool field = true,
    int uniqresid = 0) {

    field_.push_back(field);
    num_.push_back(num);
    name_.push_back(name_array);
    alterloc_.push_back(alterloc);
    resname_.push_back(resname_array);
    chain_.push_back(chain_array);
    insertres_.push_back(insertres);
    resid_.push_back(res_id);
    x_.push_back(x);
    y_.push_back(y);
    z_.push_back(z);
    occ_.push_back(occ);
    beta_.push_back(beta);
    elem_.push_back(elem);
    uniqresid_.push_back(uniqresid);

    return true;
}

vector<bool> Model::simple_select_atoms(const string &column, const vector<string> &values, const string &op) const{
    return simple_select_atoms_model(*this, column, values, op);
}

vector<bool> Model::select_tokens(const Token &tokens) const {

    vector<bool> bool_list;
    vector<bool> new_bool_list;
    string logical;
    bool not_flag = false;
    float distance = 0.0;

    // Case for simple selection
    if (tokens.is_list()) {
        const auto &token_list = tokens.as_list();
        if (!token_list.empty() && token_list[0].is_string()) {
            const string &first = token_list[0].as_string();

            // Case: simple operator-based or keyword-based selection
            if (KEYWORDS.count(first)) {
                if (token_list.size() >= 3 && token_list[1].is_string() && is_operator(token_list[1].as_string())) {
                    return simple_select_atoms(first, {token_list[2].as_string()}, token_list[1].as_string());
                } else {
                    vector<string> values;
                    values.reserve(token_list.size() - 1);
                    for (size_t i = 1; i < token_list.size(); ++i) {
                        values.push_back(token_list[i].as_string());
                    }
                    return simple_select_atoms(first, values, "isin");
                }
            }
        } 
    }

    // Nested structure
    const auto &nested_tokens = tokens.as_list();
    for (size_t i = 0; i < nested_tokens.size(); ++i) {
        const Token &tok = nested_tokens[i];

        if (tok.is_string()) {
            string token_str = tok.as_string();
            if (token_str == "and" || token_str == "or") {
                logical = token_str;
                bool_list = move(new_bool_list);
                new_bool_list.clear();
                continue;
            } else if (token_str == "not") {
                not_flag = true;
                continue;
            }
            else if (token_str == "within") {
                const auto &token_list = tokens.as_list();
                distance = stof(token_list[1].as_string());
                new_bool_list = select_tokens(token_list[3]);
                return dist_under_index(*this, new_bool_list, distance);
            }
        }

        new_bool_list = move(select_tokens(tok)); // Avoid copy

        if (not_flag) {
            for (size_t j = 0; j < new_bool_list.size(); ++j) {
                new_bool_list[j] = !new_bool_list[j];
            }
            not_flag = false;
        }

        if (!bool_list.empty() && !logical.empty()) {
            const size_t n = bool_list.size();
            if (new_bool_list.size() != n) {
                throw runtime_error("Mismatched boolean vector sizes in logical operation.");
            }
            if (logical == "and") {
                // transform(bool_list.begin(), bool_list.end(), new_bool_list.begin(), new_bool_list.begin(), logical_and<bool>());
                for (size_t k = 0; k < bool_list.size(); ++k) {
                    new_bool_list[k] = bool_list[k] && new_bool_list[k];
                }
            }
            else if (logical == "or") {
                // transform(bool_list.begin(), bool_list.end(), new_bool_list.begin(), new_bool_list.begin(), logical_or<bool>());
                for (size_t k = 0; k < bool_list.size(); ++k) {
                    new_bool_list[k] = bool_list[k] || new_bool_list[k];
                }
            }
            logical.clear();
        }
    }

    return new_bool_list;
}

vector<bool> Model::select_atoms(const string selection) const{
    Token parsed_selection = parse_selection(selection);
    return select_tokens(parsed_selection);
}

Model Model::select_index(const vector<bool> &indexes) const {
    Model selected;
    selected.clear();

    for (size_t i = 0; i < indexes.size(); ++i) {
        if (indexes[i]) {
            selected.addAtom(
                num_[i],
                name_[i],
                resname_[i],
                resid_[i],
                chain_[i],
                x_[i], y_[i], z_[i], occ_[i], beta_[i],
                alterloc_[i],
                elem_[i],
                insertres_[i],
                field_[i],
                uniqresid_[i]);
        }
    }
    return selected;
}

vector<int> Model::get_index_select(const string selection) const{
    vector<bool> bool_indexes = select_atoms(selection);
    vector<int> indices;
    for (size_t i = 0; i < bool_indexes.size(); ++i) {
        if (bool_indexes[i]) {
            indices.push_back(i);
        }
    }
    return indices;
}

vector<array<char, 2>> Model::get_uniq_chain() const{
    unordered_set<string> uniq_chains;
    vector<array<char, 2>> uniq_chain;

    for (const auto &chain : chain_) {
        string chain_str(chain.data(), strnlen(chain.data(), 2));
        if (uniq_chains.insert(chain_str).second) {
            uniq_chain.push_back(chain);
        }
    }

    return uniq_chain;

}

float Model::distance(size_t i, size_t j) const{
    return calculate_distance(x_[i], y_[i], z_[i], x_[j], y_[j], z_[j]);
}

std::array<float, 3> Model::get_centroid() const {
    if (x_.empty()) {
        throw std::runtime_error("Cannot calculate centroid of empty model");
    }
    
    float sum_x = 0.0f, sum_y = 0.0f, sum_z = 0.0f;
    
    for (size_t i = 0; i < x_.size(); ++i) {
        sum_x += x_[i];
        sum_y += y_[i];
        sum_z += z_[i];
    }
    
    float n = static_cast<float>(x_.size());
    return {sum_x / n, sum_y / n, sum_z / n};
}

std::array<float, 3> Model::get_centroid(const std::vector<int>& indices) const {
    if (indices.empty()) {
        throw std::runtime_error("Cannot calculate centroid of empty selection");
    }
    
    float sum_x = 0.0f, sum_y = 0.0f, sum_z = 0.0f;
    
    for (int idx : indices) {
        if (idx < 0 || static_cast<size_t>(idx) >= x_.size()) {
            throw std::runtime_error("Index out of bounds in centroid calculation");
        }
        sum_x += x_[idx];
        sum_y += y_[idx];
        sum_z += z_[idx];
    }
    
    float n = static_cast<float>(indices.size());
    return {sum_x / n, sum_y / n, sum_z / n};
}