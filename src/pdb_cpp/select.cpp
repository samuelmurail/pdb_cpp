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
#include <unordered_set>

#include "select.h"
#include "Model.h"
#include "geom.h"

using namespace std;


const unordered_set<string> KEYWORDS = {
    "resname", "chain", "name", "altloc",
    "resid", "residue", "beta", "occupancy",
    "x", "y", "z"
};

// NICKNAMES
unordered_map<string, string> NICKNAMES = {
    {"protein", "resname ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL"},
    {"dna", "resname DA DC DG DT"},
    {"backbone", "resname ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL and name N CA C O"},
    {"noh", "not name H*"},
    {"ions", "resname NA CL CA MG ZN MN FE CU CO NI CD K"}};


bool is_operator(const string &token)
{
    // List of valid operators
    static const vector<string> operators = {"==", "!=", ">", ">=", "<", "<="};
    // Check if the token is in the list of operators
    return find(operators.begin(), operators.end(), token) != operators.end();
}

void print_tokens(const Token &token, int indent) {
    if (token.is_string()) {
        cout << string(indent, ' ') << token.as_string() << '\n';
    } else {
        cout << string(indent, ' ') << "(\n";
        for (const auto &sub : token.as_list()) {
            print_tokens(sub, indent + 2);
        }
        cout << string(indent, ' ') << ")\n";
    }
}

// Function to replace all occurrences of a substring in a string
void replace_all(string &str, const string &from, const string &to) {
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != string::npos) {
        if (str[start_pos + 1] == '=' && (from == ">" || from =="<")) { // Handle special case for >= and <=
            start_pos += 2;
            continue;
        }
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
}

TokenList split(const string &str)
{
    istringstream iss(str);
    TokenList tokens;
    string token;
    while (iss >> token)
    {
        tokens.emplace_back(token); // Each token is stored as a string in the Token variant
    }
    return tokens;
}

// Recursive parser
TokenList parse_parentheses(const TokenList &tokens, size_t start) {
    TokenList result;
    size_t i = start;

    while (i < tokens.size()) {
        if (tokens[i].is_string() && tokens[i].as_string() == "(") {
            TokenList nested = parse_parentheses(tokens, i + 1);
            result.emplace_back(nested);
            while (i < tokens.size() && !(tokens[i].is_string() && tokens[i].as_string() == ")")) {
                ++i;
            }
        } else if (tokens[i].is_string() && tokens[i].as_string() == ")") {
            return result;
        } else {
            result.emplace_back(tokens[i]);
        }
        ++i;
    }

    return result;
}

TokenList parse_keywords(const TokenList &tokens)
{
    TokenList new_tokens;
    TokenList local_sel;

    for (const auto &token : tokens) {
        if (token.is_string()) {
            const string &val = token.as_string();

            if (KEYWORDS.count(val)) {
                if (!local_sel.empty() && local_sel[0].is_string() && local_sel[0].as_string() == "within") {
                    new_tokens.push_back(Token(local_sel));
                } else if (!local_sel.empty()) {
                    throw runtime_error("Invalid selection string: " + val);
                }
                local_sel = {token};
            } else if (val == "and" || val == "or" || val == "not") {
                if (!local_sel.empty()) {
                    new_tokens.push_back(Token(local_sel));
                    local_sel.clear();
                }
                new_tokens.push_back(token);
            } else {
                local_sel.push_back(token);
            }
        } else if (token.is_list()) {
            if (!local_sel.empty()) {
                new_tokens.push_back(Token(local_sel));
                local_sel.clear();
            }
            TokenList nested = parse_keywords(token.as_list());
            new_tokens.push_back(Token(nested));
        }
    }

    if (!local_sel.empty()) {
        new_tokens.push_back(Token(local_sel));
    }

    return new_tokens;
}


TokenList parse_within(const TokenList& tokens) {
    TokenList new_tokens;
    size_t i = 0;
    while (i < tokens.size()) {
        // if (!tokens[i].is_list()){
        //     cout << "parse_within: string : tokens[" << i << "] = " << tokens[i].as_string() << endl;
        // } else {
        //     cout << "parse_within: list :   tokens[" << i << "] = (";
        //     for (const auto& sub : tokens[i].as_list()) {
        //         cout << sub.as_string() << " ";
        //     }
        //     cout << ")" << endl;
        // }
        if (tokens[i].is_list()) {
            const auto& sublist = tokens[i].as_list();
            if (!sublist.empty() && sublist[0].is_string() && sublist[0].as_string() == "within") {
                TokenList new_token = sublist;
                if (tokens[i + 1].is_list() || (tokens[i + 1].is_string() && tokens[i + 1].as_string() != "not")) {
                    new_token.emplace_back(tokens[i + 1]);
                    ++i;
                } else {
                    new_token.emplace_back(Token(TokenList{tokens[i + 1], tokens[i + 2]}));
                    i += 2;
                }
                new_tokens.push_back(Token(new_token));
            } else {
                new_tokens.push_back(Token(parse_within(sublist)));
            }
        } else {
            new_tokens.push_back(tokens[i]);
        }
        ++i;
    }
    return new_tokens;
}


// Main function to parse the selection string
Token parse_selection(string selection)
{

    // Add spaces around operators and parentheses
    for (const string &op : initializer_list<string>{"!=", "==", "<=", ">=", ":", "and", "or", "not", "(", ")", "<", ">", }) {
        replace_all(selection, op, " " + op + " ");
    }

    // Replace nicknames with their corresponding values
    for (const auto &pair : NICKNAMES) {
        replace_all(selection, pair.first, pair.second);
    }

    TokenList tokens = split(selection);
    // cout << "tokens before parentheses: " << endl;
    // print_tokens(tokens);
    tokens = parse_parentheses(tokens);
    // cout << "tokens after parentheses: " << endl;
    // print_tokens(tokens);
    tokens = parse_keywords(tokens);
    // cout << "tokens after keywords: " << endl;
    // print_tokens(tokens);
    // Parse within keyword
    // cout << "tokens after keywords: " << endl;
    // print_tokens(tokens);
    tokens = parse_within(tokens);
    // cout << "tokens after within: " << endl;
    // print_tokens(tokens);

    return tokens;
}

vector<bool> simple_select_atoms_model(const Model &model, const string &column, const vector<string> &values, const string &op)
{
    size_t n = model.size();
    vector<bool> result(n, false);

    if (0) {
        cout << "column: " << column << endl;
        cout << "op: " << op << endl;
        cout << "values: ";
        for (const auto &v : values) {
            cout << v << " ";
        }
        cout << endl;
    }
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
        const vector<float> &model_val = (column == "x") ? model.get_x() : (column == "y") ? model.get_y()
                                                                       : (column == "z")   ? model.get_z()
                                                                       : (column == "occ") ? model.get_occ()
                                                                                           : model.get_beta();
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
            for (const auto &v : values)
                value_set.insert(stof(v));
            for (size_t i = 0; i < n; ++i)
                result[i] = value_set.count(model_val[i]) > 0;
            return result;
        } else {
            throw invalid_argument("Unsupported operator for '" + column + "'");
        }
    } else if (column == "resid" || column == "residue" || column == "num") {
        const vector<int> &model_val = (column == "resid") ? model.get_resid() : (column == "residue") ? model.get_uniqresid()
                                                                                                         : model.get_num();
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
        } else if (op == "isin") {
            set<int> value_set;
            for (const auto &v : values)
                value_set.insert(stoi(v));
            for (size_t i = 0; i < n; ++i)
                result[i] = value_set.count(model_val[i]) > 0;
            return result;
        } else {
            throw invalid_argument("Unsupported operator for '" + column + "'");
        }
    } else if (column == "chain" || column == "altloc" || column == "insertres") {
        const vector<array<char, 2>> &model_val = (column == "chain") ? model.get_chain() : (column == "altloc") ? model.get_alterloc()
                                                                                                                 : model.get_insertres();
        if (op == "==") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i][0] == values[0][0];
            }
        } else if (op == "!=") {
            for (size_t i = 0; i < n; ++i) {
                result[i] = model_val[i][0] != values[0][0];
            }
        }
        else if (op == "isin") {
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

vector<bool> dist_under_index(Model &model, vector<bool> selection, float distance) {
    vector<bool> result(selection.size(), false);
    vector<float> x = model.get_x();
    vector<float> y = model.get_y();
    vector<float> z = model.get_z();
    float square_distance = distance * distance, sq_dist;

    for (size_t i = 0; i < selection.size(); ++i) {
        if (selection[i]) {
            result[i] = true;
        } else {
            for (size_t j = 0; j < selection.size(); ++j) {
                if (selection[j]) {
                    sq_dist = calculate_square_distance(x[i], y[i], z[i], x[j], y[j], z[j]);
                    if (sq_dist < square_distance) {
                        result[i] = true;
                        break;
                    }
                }
            }
        }
    }
    return result;
}