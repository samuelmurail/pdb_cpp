#include <array>
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

    // If not found, throw an exception or return a placeholder
    throw invalid_argument("Unknown residue name: " + resname);
}