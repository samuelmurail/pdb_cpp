#include <cstring>
#include <iomanip>
#include <set>

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
    size_t n = this->size();
    vector<bool> result(n, false);

    cout << "column: " << column << endl;
    cout << "op: " << op << endl;
    cout << "values: ";
    for (const auto &v : values) {
        cout << v << " ";
    }

    auto str_equal = [](const array<char, 5> &a, const string &b) {
        return strncmp(a.data(), b.c_str(), 5) == 0;
    };

    auto str_startswith = [](const array<char, 5> &a, const string &b) {
        return strncmp(a.data(), b.c_str(), b.size()) == 0;
    };

    if (column == "name") {
        if (op == "==") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = str_equal(name_[i], values[0]);
            }
        } else if (op == "!=") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = !str_equal(name_[i], values[0]);
            }
        } else if (op == "startswith") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = str_startswith(name_[i], values[0]);
            }
        } else if (op == "isin") {
            set<string> value_set(values.begin(), values.end());
            for (size_t i = 0; i < n; ++i) {
                string val(name_[i].data(), strnlen(name_[i].data(), 5));
                result[i] = value_set.count(val) > 0;
            }
        } else {
            throw invalid_argument("Unsupported operator for 'name'");
        }

    } else if (column == "resname") {
        if (op == "==") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = str_equal(resname_[i], values[0]);
            }
        } else if (op == "isin") {
            set<string> value_set(values.begin(), values.end());
            for (size_t i = 0; i < n; ++i) {
                string val(resname_[i].data(), strnlen(resname_[i].data(), 5));
                result[i] = value_set.count(val) > 0;
            }
        } else {
            throw invalid_argument("Unsupported operator for 'resname'");
        }

    } else if (column == "chain") {
        if (op == "==") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = chain_[i][0] == values[0][0];
            }
        } else if (op == "!=") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = chain_[i][0] != values[0][0];
            }
        } else {
            throw invalid_argument("Unsupported operator for 'chain'");
        }

    } else if (column == "resid") {
        int val = stoi(values[0]);
        if (op == "==") {
            for (size_t i = 0; i < n; ++i) result[i] = resid_[i] == val;
        } else if (op == "!=") {
            for (size_t i = 0; i < n; ++i) result[i] = resid_[i] != val;
        } else if (op == "<") {
            for (size_t i = 0; i < n; ++i) result[i] = resid_[i] < val;
        } else if (op == "<=") {
            for (size_t i = 0; i < n; ++i) result[i] = resid_[i] <= val;
        } else if (op == ">") {
            for (size_t i = 0; i < n; ++i) result[i] = resid_[i] > val;
        } else if (op == ">=") {
            for (size_t i = 0; i < n; ++i) result[i] = resid_[i] >= val;
        } else if (op == "isin") {
            set<int> valset;
            for (const auto &v : values) valset.insert(stoi(v));
            for (size_t i = 0; i < n; ++i) result[i] = valset.count(resid_[i]) > 0;
        } else {
            throw invalid_argument("Unsupported operator for 'resid'");
        }

    } else if (column == "x") {
        float val = stof(values[0]);
        const vector<float> &data = x_;
        if (op == "==") {
            for (size_t i = 0; i < n; ++i) result[i] = data[i] == val;
        } else if (op == "!=") {        } else if (op == "!=") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = !str_equal(name_[i], values[0]);
            }
        } else if (op == "startswith") {

            for (size_t i = 0; i < n; ++i) result[i] = data[i] != val;
        } else if (op == "<") {
            for (size_t i = 0; i < n; ++i) result[i] = data[i] < val;
        } else if (op == "<=") {
            for (size_t i = 0; i < n; ++i) result[i] = data[i] <= val;
        } else if (op == ">") {
            for (size_t i = 0; i < n; ++i) result[i] = data[i] > val;
        } else if (op == ">=") {
            for (size_t i = 0; i < n; ++i) result[i] = data[i] >= val;
        } else {
            throw invalid_argument("Unsupported operator for 'x'");
        }

    } else {
        throw invalid_argument("Column '" + column + "' not recognized");
    }

    return result;
}


// vector<bool> Model::simple_select_atoms(const string &column, const vector<string> &values, const string &op = "==") {
//     size_t n = this->size();
//     vector<bool> result(n, false); // Initialize a boolean vector with the size of the model

//     // Lambda function to compare based on the operator
//     auto compare = [&](const auto &atom_value, const string &value) -> bool {
//         if (op == "==") return atom_value == value;
//         if (op == "!=") return atom_value != value;
//         if (op == ">") return stof(atom_value) > stof(value);
//         if (op == ">=") return stof(atom_value) >= stof(value);
//         if (op == "<") return stof(atom_value) < stof(value);
//         if (op == "<=") return stof(atom_value) <= stof(value);
//         if (op == "isin") return find(values.begin(), values.end(), atom_value) != values.end();
//         throw invalid_argument("Invalid operator: " + op);
//     };

// // Select based on the column
// if (column == "name") {
//     for (size_t i = 0; i < name_.size(); ++i) {
//         for (const auto &value : values) {
//             if (compare(name_[i], value)) {
//                 result[i] = true;
//                 break;
//             }
//         }
//     }
// }
// } else if (column == "chain") {
//     for (size_t i = 0; i < chain_.size(); ++i) {
//         for (const auto &value : values) {
//             if (compare(chain_[i], value)) {
//                 result[i] = true;
//                 break;
//             }
//         }
//     }
// } else if (column == "name") {
//     for (size_t i = 0; i < name_.size(); ++i) {
//         for (const auto &value : values) {
//             if (compare(name_[i], value)) {
//                 result[i] = true;
//                 break;
//             }
//         }
//     }
// } else if (column == "altloc") {
//     for (size_t i = 0; i < alterloc_.size(); ++i) {
//         for (const auto &value : values) {
//             if (compare(alterloc_[i], value)) {
//                 result[i] = true;
//                 break;
//             }
//         }
//     }
// } else if (column == "resid") {
//     for (size_t i = 0; i < resid_.size(); ++i) {
//         for (const auto &value : values) {
//             if (compare(to_string(resid_[i]), value)) {
//                 result[i] = true;
//                 break;
//             }
//         }
//     }
// } else if (column == "beta") {
//     for (size_t i = 0; i < beta_.size(); ++i) {
//         for (const auto &value : values) {
//             if (compare(to_string(beta_[i]), value)) {
//                 result[i] = true;
//                 break;
//             }
//         }
//     }
// } else if (column == "occ") {
//     for (size_t i = 0; i < occ_.size(); ++i) {
//         for (const auto &value : values) {
//             if (compare(to_string(occ_[i]), value)) {
//                 result[i] = true;
//                 break;
//             }
//         }
//     }
// } else if (column == "x" || column == "y" || column == "z") {
//     const auto &coords = (column == "x") ? x_ : (column == "y") ? y_ : z_;
//     for (size_t i = 0; i < coords.size(); ++i) {
//         for (const auto &value : values) {
//             if (compare(to_string(coords[i]), value)) {
//                 result[i] = true;
//                 break;
//             }
//         }
//     }
// } else {
//     throw invalid_argument("Invalid column: " + column);
// }
//     return result;
// }


vector<bool> Model::select_tokens(const Token &tokens) {
    vector<bool> bool_list;
    vector<bool> new_bool_list;
    string logical;
    bool not_flag = false;

    // // Case for simple selection
    // if (is_simple_list(tokens)) {
    //     const auto &token_list = get<vector<Token>>(tokens.value);
    //     if (token_list.size() >= 3 && is_operator(get<string>(token_list[1].value))) {
    //         return simple_select_atoms(
    //             get<string>(token_list[0].value), // column
    //             get<string>(token_list[2].value), // values
    //             get<string>(token_list[1].value)  // operator
    //         );
    //     } else {
    //         vector<string> values;
    //         for (size_t i = 1; i < token_list.size(); ++i) {
    //             values.push_back(get<string>(token_list[i].value));
    //         }
    //         return simple_select_atoms(get<string>(token_list[0].value), values);
    //     }
    // }
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

    // // Process nested tokens
    // const auto &token_list = get<vector<Token>>(tokens.value);
    // for (size_t i = 0; i < token_list.size(); ++i) {
    //     if (holds_alternative<string>(token_list[i].value)) {
    //         string token = get<string>(token_list[i].value);

    //         if (token == "and" || token == "or") {
    //             logical = token;
    //             bool_list = new_bool_list;
    //             new_bool_list.clear();
    //             continue;
    //         } else if (token == "not") {
    //             not_flag = true;
    //             continue;
    //         }
    //     }

    //     // Recursive call for nested tokens
    //     new_bool_list = select_tokens(token_list[i]);

    //     if (not_flag) {
    //         for (auto &val : new_bool_list) {
    //             val = !val;
    //         }
    //         not_flag = false;
    //     }

    //     if (!bool_list.empty() && !logical.empty()) {
    //         if (logical == "and") {
    //             transform(bool_list.begin(), bool_list.end(), new_bool_list.begin(), new_bool_list.begin(), logical_and<bool>());
    //         } else if (logical == "or") {
    //             transform(bool_list.begin(), bool_list.end(), new_bool_list.begin(), new_bool_list.begin(), logical_or<bool>());
    //         }
    //         logical.clear();
    //     }
    // }

    return new_bool_list;
}