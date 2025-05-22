#ifndef SELECT_H
#define SELECT_H

#include <vector>
#include <string>
#include <unordered_map>
#include <variant>

using namespace std;

class Model;

// Define a recursive variant: a token can be either a string or a nested list of tokens
struct Token;
using TokenVariant = variant<string, vector<Token>>;
struct Token {
    TokenVariant value;

    // Constructors for convenience
    Token() = default;
    Token(const string& str) : value(str) {}
    Token(const vector<Token>& tokens) : value(tokens) {}
};

// Alias for list of tokens
using TokenList = vector<Token>;

// Function declarations
bool is_simple_list(const Token &tokens);
bool is_operator(const string &token);
void print_tokens(const Token& token, int indent=0);
void replace_all(string &str, const string &from, const string &to);
TokenList split(const string &str);
TokenList parse_parentheses(const TokenList &tokens, size_t start);
Token parse_selection(string selection, const unordered_map<string, string> &nicknames);

vector<bool> simple_select_atoms_model(const Model &model, const string &column, const vector<string> &values, const string &op);

// Keywords and nicknames
extern std::vector<std::string> KEYWORDS;
extern unordered_map<string, string> NICKNAMES;

#endif // SELECT_H