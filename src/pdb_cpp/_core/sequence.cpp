#include <array>
#include <cctype>
#include <cstring>
#include <stdexcept>

#include "data/residue.h"

using namespace std;

char convert_to_one_letter_resname(const array<char, 5> &resname_array) {
    // Convert the array to a string (trim trailing null characters)
    string resname(resname_array.data(), strnlen(resname_array.data(), resname_array.size()));

    // Look up the residue name in the dictionary
    auto it = AA_DICT_L.find(resname);
    if (it != AA_DICT_L.end()) {
        return it->second; // Return the 1-letter code
    }

    // // If not found, throw an exception or return a placeholder
    // throw invalid_argument("Unknown residue name: " + resname);

    // If not found, return a placeholder (e.g., 'X' for unknown)
    return 'X';
}

char convert_to_one_letter_resname_dl(const array<char, 5> &resname_array) {
    string resname(resname_array.data(), strnlen(resname_array.data(), resname_array.size()));

    auto it_l = AA_DICT_L.find(resname);
    if (it_l != AA_DICT_L.end()) {
        return it_l->second;
    }

    auto it_d = AA_DICT_D.find(resname);
    if (it_d != AA_DICT_D.end()) {
        return static_cast<char>(tolower(static_cast<unsigned char>(it_d->second)));
    }

    auto it_na = NA_DICT.find(resname);
    if (it_na != NA_DICT.end()) {
        return it_na->second;
    }

    // throw invalid_argument("Unknown residue name: " + resname);

    // If not found, return a placeholder (e.g., 'X' for unknown)
    return 'X';
}

char convert_to_one_letter_resname_na(const array<char, 5> &resname_array) {
    string resname(resname_array.data(), strnlen(resname_array.data(), resname_array.size()));

    auto it = AA_NA_DICT.find(resname);
    if (it != AA_NA_DICT.end()) {
        return it->second;
    }

    // If not found, return a placeholder (e.g., 'X' for unknown)
    return 'X';
}

char convert_to_one_letter_resname_any(const array<char, 5> &resname_array) {
    string resname(resname_array.data(), strnlen(resname_array.data(), resname_array.size()));

    auto it_l = AA_DICT_L.find(resname);
    if (it_l != AA_DICT_L.end()) {
        return it_l->second;
    }

    auto it_d = AA_DICT_D.find(resname);
    if (it_d != AA_DICT_D.end()) {
        return it_d->second;
    }

    //throw invalid_argument("Unknown residue name: " + resname);
    // If not found, return a placeholder (e.g., 'X' for unknown)
    return 'X';
}