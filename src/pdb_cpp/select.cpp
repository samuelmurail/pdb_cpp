#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <stdexcept>
#include <iostream>
#include <variant>
#include <set>

#include "select.h"
#include "Model.h"


using namespace std;

std::vector<std::string> KEYWORDS = {
    "resname",
    "chain",
    "name",
    "altloc",
    "resid",
    "residue",
    "beta",
    "occ",
    "x",
    "y",
    "z"
};

//NICKNAMES
unordered_map<string, string> NICKNAMES = {
    {"protein", "resname ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL"},
    {"dna", "resname DA DC DG DT"},
    {"backbone", "resname ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL and name N CA C O"},
    {"noh", "not name H*"},
    {"ions", "resname NA CL CA MG ZN MN FE CU CO NI CD K"}
};

bool is_simple_list(const Token &tokens) {
    // Check if the token is a list of tokens
    if (!holds_alternative<vector<Token>>(tokens.value)) {
        return false;
    }

    // Get the list of tokens
    const auto &token_list = get<vector<Token>>(tokens.value);

    // Check if any token in the list is itself a nested list
    for (const auto &token : token_list) {
        if (holds_alternative<vector<Token>>(token.value)) {
            return false;
        }
    }
    return true;
}

bool is_operator(const string &token) {
    // List of valid operators
    static const vector<string> operators = {"==", "!=", ">", ">=", "<", "<="};
    // Check if the token is in the list of operators
    return find(operators.begin(), operators.end(), token) != operators.end();
}

void print_tokens(const Token& token, int indent) {
    if (holds_alternative<string>(token.value)) {
        cout << string(indent, ' ') << get<string>(token.value) << '\n';
    } else {
        cout << string(indent, ' ') << "(\n";
        for (const auto& sub : get<vector<Token>>(token.value)) {
            print_tokens(sub, indent + 2);
        }
        cout << string(indent, ' ') << ")\n";
    }
}

// Function to replace all occurrences of a substring in a string
void replace_all(string &str, const string &from, const string &to) {
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
}

TokenList split(const string &str) {
    istringstream iss(str);
    TokenList tokens;
    string token;
    while (iss >> token) {
        tokens.emplace_back(token); // Each token is stored as a string in the Token variant
    }
    return tokens;
}

// Recursive parser
TokenList parse_parentheses(const TokenList &tokens, size_t start = 0) {
    TokenList result;
    size_t i = start;

    while (i < tokens.size()) {
        if (holds_alternative<string>(tokens[i].value) && get<string>(tokens[i].value) == "(") {
            // Recursively parse the nested parentheses
            TokenList nested = parse_parentheses(tokens, i + 1);
            result.emplace_back(nested); // Add the nested TokenList as a Token
            // Skip to the position after the closing parenthesis
            while (i < tokens.size() && !(holds_alternative<string>(tokens[i].value) && get<string>(tokens[i].value) == ")")) {
                ++i;
            }
        } else if (holds_alternative<string>(tokens[i].value) && get<string>(tokens[i].value) == ")") {
            // End of the current nested list
            return result;
        } else {
            // Add a simple string token
            result.emplace_back(tokens[i]);
        }
        ++i;
    }

    return result;
}

// Main function to parse the selection string
Token parse_selection(string selection, const unordered_map<string, string> &nicknames) {

    // Add spaces around operators and parentheses
    for (const string &op : initializer_list<string>{"(", ")", "<", ">", "!=", "==", "<=", ">=", ":", "and", "or", "not"}) {
        replace_all(selection, op, " " + op + " ");
    }

    // Replace nicknames with their corresponding values
    for (const auto &pair : nicknames) {
        replace_all(selection, pair.first, pair.second);
    }

    TokenList tokens = split(selection);

    tokens = parse_parentheses(tokens);

    return tokens;
}

vector<bool> simple_select_atoms_model(const Model &model, const string &column, const vector<string> &values, const string &op) {
    size_t n = model.size();
    vector<bool> result(n, false);

    cout << "column: " << column << endl;
    cout << "op: " << op << endl;
    cout << "values: ";
    for (const auto &v : values) {
        cout << v << " ";
    }
    cout << endl;

    auto str_equal = [](const array<char, 5> &a, const string &b) {
        return strncmp(a.data(), b.c_str(), 5) == 0;
    };

    auto str_startswith = [](const array<char, 5> &a, const string &b) {
        return strncmp(a.data(), b.c_str(), b.size()) == 0;
    };

    if (column == "name" || column == "resname" || column == "elem") {
        vector<array<char, 5>> model_val = (column == "name") ? model.get_name() : (column == "resname") ? model.get_resname() : model.get_elem();

        if (op == "==") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = str_equal(model_val[i], values[0]);
            }
            return result;
        } else if (op == "!=") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = !str_equal(model_val[i], values[0]);
            }
            return result;
        } else if (op == "startswith") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = str_startswith(model_val[i], values[0]);
            }
            return result;
        } else if (op == "isin") {
            set<string> value_set(values.begin(), values.end());
            for (size_t i = 0; i < n; ++i) {
                string val(model_val[i].data(), strnlen(model_val[i].data(), 5));
                result[i] = value_set.count(val) > 0;
            }
            return result;
        } else {
            throw invalid_argument("Unsupported operator for 'name'");
        }
    } else if (column == "x" || column == "y" || column == "z" || column == "occ" || column == "beta") {
        const vector<float> &model_val = (column == "x") ? model.get_x() : (column == "y") ? model.get_y() : (column == "z") ? model.get_z() : (column == "occ") ? model.get_occ() : model.get_beta();
        float val = stof(values[0]);
        if (op == "==") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i] == val;
            }
            return result;
        } else if (op == "!=") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i] != val;
            }
            return result;
        } else if (op == "<") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i] < val;
            }
            return result;
        } else if (op == "<=") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i] <= val;
            }
            return result;
        } else if (op == ">") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i] > val;
            }
            return result;
        } else if (op == ">=") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i] >= val;
            }
            return result;
        } else if (op == "isin") {
            set<float> value_set;
            for (const auto &v : values) value_set.insert(stof(v));
            for (size_t i = 0; i < n; ++i) result[i] = value_set.count(model_val[i]) > 0;
            return result;
        } else {
            throw invalid_argument("Unsupported operator for '" + column + "'");
        }
    } else if (column == "resid" || column == "uniqresid" || column == "num") {
        const vector<int> &model_val = (column == "resid") ? model.get_resid() : (column == "uniqresid") ? model.get_uniqresid() : model.get_num();
        int val = stoi(values[0]);
        if (op == "==") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i] == val;
            }
            return result;
        } else if (op == "!=") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i] != val;
            }
            return result;
        } else if (op == "<") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i] < val;
            }
            return result;
        } else if (op == "<=") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i] <= val;
            }
            return result;
        } else if (op == ">") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i] > val;
            }
            return result;
        } else if (op == ">=") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i] >= val;
            }
            return result;
        } else
        if (op == "isin") {
            set<int> value_set;
            for (const auto &v : values) value_set.insert(stoi(v));
            for (size_t i = 0; i < n; ++i) result[i] = value_set.count(model_val[i]) > 0;
            return result;
        } else {
            throw invalid_argument("Unsupported operator for '" + column + "'");
        }
    } else if (column == "chain" || column == "altloc" || column == "insertres") {
        const vector<array<char, 2>> &model_val = (column == "chain") ? model.get_chain() : (column == "altloc") ? model.get_alterloc() : model.get_insertres();
        if (op == "==") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i][0] == values[0][0];
            }
        } else if (op == "!=") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i][0] != values[0][0];
            }
        } else if (op == "isin") {
            set<string> value_set(values.begin(), values.end());
            for (size_t i = 0; i < n; ++i) {
                string val(model_val[i].data(), strnlen(model_val[i].data(), 2));
                result[i] = value_set.count(val) > 0;
            }
        } else {
            throw invalid_argument("Unsupported operator for '" + column + "'");
        }
    } else {
        throw invalid_argument("Invalid column: " + column);
    }

    return result;
}


// int main() {
//     //string selection = "protein and resname ALA and chain A and (x > 10 or y < 5) and not (z == 0)";
//     string selection = "resname ALA and resid 4:10 and chain A B and (x > 10 or y < 5 and (z<10 and z>-10)) and not (z == 0)";
//     Token parsed_selection = parse_selection(selection, NICKNAMES);

//     cout << "Parsed selection: ";

//     print_tokens(parsed_selection);


//     return 0;
// }