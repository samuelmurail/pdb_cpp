#ifndef SELECT_H
#define SELECT_H

#include <vector>
#include <string>
#include <unordered_map>
#include <variant>
#include <unordered_set>

using namespace std;

class Model;


// Recursive token structure
struct Token;
using TokenVariant = variant<string, vector<Token>>;
struct Token {
    TokenVariant value;

    Token() = default;
    Token(const string& str) : value(str) {}
    Token(const vector<Token>& tokens) : value(tokens) {}

    bool is_string() const {
        return holds_alternative<string>(value);
    }
    bool is_list() const {
        return holds_alternative<vector<Token>>(value);
    }

    const string& as_string() const {
        return get<string>(value);
    }
    const vector<Token>& as_list() const {
        return get<vector<Token>>(value);
    }

    const string get_first() const {
        if (is_list() && !as_list().empty()) {
            return as_list()[0].as_string();
        }
        return "";
    }
};

// Alias for list of tokens
using TokenList = vector<Token>;


// Function declarations
bool is_simple_list(const Token &tokens);
bool is_operator(const string &token);
void print_tokens(const Token& token, int indent=0);
void replace_all(string &str, const string &from, const string &to);
TokenList split(const string &str);
TokenList parse_parentheses(const TokenList &tokens, size_t start = 0);
Token parse_selection(string selection);

vector<bool> simple_select_atoms_model(const Model &model, const string &column, const vector<string> &values, const string &op);
vector<bool> dist_under_index(const Model &model, vector<bool> selection, float distance);

// Keywords and nicknames
extern const unordered_set<string> KEYWORDS;
extern unordered_map<string, string> NICKNAMES;

#endif // SELECT_H