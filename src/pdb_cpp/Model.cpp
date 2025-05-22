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

    cout << "Selecting tokens..." << endl;


    // Case for simple selection
    if (is_simple_list(tokens)) {
        const auto &token_list = get<vector<Token>>(tokens.value);
        if (token_list.size() >= 3 && is_operator(get<string>(token_list[1].value))) {
            cout << "Simple selection with operator: " << get<string>(token_list[1].value) << endl;
            // values for operators like ==, !=, etc. should be a vector, not a single string
            std::vector<std::string> values = {get<string>(token_list[2].value)};
            return simple_select_atoms(
                get<string>(token_list[0].value), // column
                values, // values as vector
                get<string>(token_list[1].value)  // operator
            );
        } else {
            cout << "Simple selection without operator: " << get<string>(token_list[0].value) << endl;
            vector<string> values;
            for (size_t i = 1; i < token_list.size(); ++i) {
                values.push_back(get<string>(token_list[i].value));
            }
            // Provide default operator for simple selection (e.g., "isin")
            return simple_select_atoms(get<string>(token_list[0].value), values, "isin");
        }
    }
    // // Case for "within" selection
    // else if (holds_alternative<string>(tokens.value) && get<string>(tokens.value) == "within") {
    //     const auto &token_list = get<vector<Token>>(tokens.value);
    //     if (token_list.size() != 4) {
    //         throw invalid_argument("within selection must have 3 arguments");
    //     }
    //     new_bool_list = select_tokens(token_list[3]);
    //     float distance = stof(get<string>(token_list[1].value));
    //     vector<int> sel_2 = select_index(new_bool_list);

    //     return dist_under_index(sel_2, distance);
    // }

    // Process nested tokens
    const auto &token_list_nested = get<vector<Token>>(tokens.value);
    for (size_t i = 0; i < token_list_nested.size(); ++i) {
        cout << "Processing token: " << i << endl;
        cout << "Token value: " << get<string>(token_list_nested[i].value) << endl;
        if (holds_alternative<string>(token_list_nested[i].value)) {
            string token = get<string>(token_list_nested[i].value);

            if (token == "and" || token == "or") {
                logical = token;
                bool_list = new_bool_list;
                new_bool_list.clear();
                continue;
            } else if (token == "not") {
                not_flag = true;
                continue;
            }
        }

        // Recursive call for nested tokens
        new_bool_list = select_tokens(token_list_nested[i]);

        if (not_flag) {
            for (size_t j = 0; j < new_bool_list.size(); ++j) {
                new_bool_list[j] = !new_bool_list[j];
            }
            not_flag = false;
        }

        if (!bool_list.empty() && !logical.empty()) {
            if (logical == "and") {
                transform(bool_list.begin(), bool_list.end(), new_bool_list.begin(), new_bool_list.begin(), logical_and<bool>());
            } else if (logical == "or") {
                transform(bool_list.begin(), bool_list.end(), new_bool_list.begin(), new_bool_list.begin(), logical_or<bool>());
            }
            logical.clear();
        }
    }

    return new_bool_list;
}