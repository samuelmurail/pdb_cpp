#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <stdexcept>
#include <iostream>
#include <variant>

#include "select.h"

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

void print_tokens(const Token& token, int indent = 0) {
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

// int main() {
//     //string selection = "protein and resname ALA and chain A and (x > 10 or y < 5) and not (z == 0)";
//     string selection = "resname ALA and resid 4:10 and chain A B and (x > 10 or y < 5 and (z<10 and z>-10)) and not (z == 0)";
//     Token parsed_selection = parse_selection(selection, NICKNAMES);

//     cout << "Parsed selection: ";

//     print_tokens(parsed_selection);


//     return 0;
// }