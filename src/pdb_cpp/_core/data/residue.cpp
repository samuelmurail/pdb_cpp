#include <unordered_map>
#include <string>

using namespace std;

// Amino acid dictionary (L)
unordered_map<string, char> AA_DICT_L = {
    {"GLY", 'G'},
    {"HIS", 'H'},
    {"HSP", 'H'},
    {"HSE", 'H'},
    {"HSD", 'H'},
    {"HIP", 'H'},
    {"HIE", 'H'},
    {"HID", 'H'},
    {"ARG", 'R'},
    {"LYS", 'K'},
    {"ASP", 'D'},
    {"ASPP",'D'},
    {"ASH", 'D'},
    {"GLU", 'E'},
    {"GLUP",'E'},
    {"GLH", 'E'},
    {"SER", 'S'},
    {"THR", 'T'},
    {"ASN", 'N'},
    {"GLN", 'Q'},
    {"CYS", 'C'},
    {"SEC", 'U'},
    {"PRO", 'P'},
    {"ALA", 'A'},
    {"ILE", 'I'},
    {"PHE", 'F'},
    {"TYR", 'Y'},
    {"TRP", 'W'},
    {"VAL", 'V'},
    {"LEU", 'L'},
    {"MET", 'M'}
};

// D-amino acid dictionary
unordered_map<string, char> AA_DICT_D = {
    {"DAL", 'A'},
    {"DAR", 'R'},
    {"DSG", 'N'},
    {"DAS", 'D'},
    {"DCY", 'C'},
    {"DGN", 'Q'},
    {"DGL", 'E'},
    {"DHI", 'H'},
    {"DIL", 'I'},
    {"DLE", 'L'},
    {"DLY", 'K'},
    {"DME", 'M'},
    {"MED", 'M'},
    {"DPH", 'F'},
    {"DPN", 'F'},
    {"DPR", 'P'},
    {"DSE", 'S'},
    {"DSN", 'S'},
    {"DTH", 'T'},
    {"DTR", 'W'},
    {"DTY", 'Y'},
    {"DVA", 'V'}
};

// Nucleic acid dictionary
unordered_map<string, char> NA_DICT = {
    {"DA", 'A'},
    {"DT", 'T'},
    {"DC", 'C'},
    {"DG", 'G'}
};

// Combined amino acid and nucleic acid dictionary
unordered_map<string, char> AA_NA_DICT = []() {
    unordered_map<string, char> combined = AA_DICT_L;
    combined.insert(AA_DICT_D.begin(), AA_DICT_D.end());
    combined.insert(NA_DICT.begin(), NA_DICT.end());
    return combined;
}();