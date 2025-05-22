#include <cstring>
#include <iomanip>

#include "Model.h"
#include "select.h"

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
           const array<char, 5>& name_array,
           const array<char, 5>& resname_array,
           int res_id,
           const array<char, 2>& chain_array,
           float x, float y, float z, float occ, float beta,
           const array<char, 2>& alterloc = {' ', '\0'},
           const array<char, 5>& elem = {' ', ' ', ' ', ' ', '\0'},
           const array<char, 2>& insertres = {' ', '\0'},
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

vector<bool> Model::simple_select_atoms(const string &column, const vector<string> &values, const string &op) {
    return simple_select_atoms_model(*this, column, values, op);
}

vector<bool> Model::select_tokens(const Token &tokens) {
    vector<bool> bool_list;
    vector<bool> new_bool_list;
    string logical;
    bool not_flag = false;

    // Case for simple selection
    if (tokens.is_list()) {
        const auto &token_list = tokens.as_list();
        if (!token_list.empty() && token_list[0].is_string()) {
            const string &first = token_list[0].as_string();

            // Case: simple operator-based or keyword-based selection
            if (KEYWORDS.count(first)) {
                if (token_list.size() >= 3 && token_list[1].is_string() && is_operator(token_list[1].as_string())) {
                    vector<string> values = {token_list[2].as_string()};
                    return simple_select_atoms(first, values, token_list[1].as_string());
                } else {
                    vector<string> values;
                    for (size_t i = 1; i < token_list.size(); ++i) {
                        values.push_back(token_list[i].as_string());
                    }
                    return simple_select_atoms(first, values, "isin");
                }
            }

            // Case: within keyword
            if (first == "within" && token_list.size() == 4) {
                float distance = stof(token_list[1].as_string());
                vector<bool> selection = select_tokens(token_list[3]);
                //vector<int> sel_index = select_index(selection);
                //return dist_under_index(sel_index, distance);
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
                bool_list = new_bool_list;
                new_bool_list.clear();
                continue;
            } else if (token_str == "not") {
                not_flag = true;
                continue;
            }
        }

        new_bool_list = select_tokens(tok);

        if (not_flag) {
            for (bool val : new_bool_list) {
                val = !val;
            }
            not_flag = false;
        }

        if (!bool_list.empty() && !logical.empty()) {
            if (logical == "and") {
                //transform(bool_list.begin(), bool_list.end(), new_bool_list.begin(), new_bool_list.begin(), logical_and<bool>());
                for (size_t k = 0; k < bool_list.size(); ++k) {
                    new_bool_list[k] = bool_list[k] && new_bool_list[k];
                }
            } else if (logical == "or") {
                transform(bool_list.begin(), bool_list.end(), new_bool_list.begin(), new_bool_list.begin(), logical_or<bool>());
            }
            logical.clear();
        }
    }

    return new_bool_list;
}

vector<bool> Model::select_atoms(const string selection) {

    Token parsed_selection = parse_selection(selection);
    return select_tokens(parsed_selection);

}