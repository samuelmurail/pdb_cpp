#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include <vector>
#include <string>
#include <array>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <memory>
#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <utility>

#include "sequence_align.h"
#include "Coor.h"

using namespace std;

#define MATRIX_SIZE 23
const std::string DEFAULT_MATRIX_FILE = "./data/blosum62.txt";

const std::unordered_map<char, std::unordered_map<char, int>> BLOSUM62 = {
    {'A', {{'A', 4}, {'R', -1}, {'N', -2}, {'D', -2}, {'C', 0}, {'Q', -1}, {'E', -1}, {'G', 0}, {'H', -2}, {'I', -1}, {'L', -1}, {'K', -1}, {'M', -1}, {'F', -2}, {'P', -1}, {'S', 1}, {'T', 0}, {'W', -3}, {'Y', -2}, {'V', 0}}},
    {'R', {{'A', -1}, {'R', 5}, {'N', 0}, {'D', -2}, {'C', -3}, {'Q', 1}, {'E', 0}, {'G', -2}, {'H', 0}, {'I', -3}, {'L', -2}, {'K', 2}, {'M', -1}, {'F', -3}, {'P', -2}, {'S', -1}, {'T', -1}, {'W', -3}, {'Y', -2}, {'V', -3}}},
    {'N', {{'A', -2}, {'R', 0}, {'N', 6}, {'D', 1}, {'C', -3}, {'Q', 0}, {'E', 0}, {'G', 0}, {'H', 1}, {'I', -3}, {'L', -3}, {'K', 0}, {'M', -2}, {'F', -3}, {'P', -2}, {'S', 1}, {'T', 0}, {'W', -4}, {'Y', -2}, {'V', -3}}},
    {'D', {{'A', -2}, {'R', -2}, {'N', 1}, {'D', 6}, {'C', -3}, {'Q', 0}, {'E', 2}, {'G', -1}, {'H', -1}, {'I', -3}, {'L', -4}, {'K', -1}, {'M', -3}, {'F', -3}, {'P', -1}, {'S', 0}, {'T', -1}, {'W', -4}, {'Y', -3}, {'V', -3}}},
    {'C', {{'A', 0}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', 9}, {'Q', -3}, {'E', -4}, {'G', -3}, {'H', -3}, {'I', -1}, {'L', -1}, {'K', -3}, {'M', -1}, {'F', -2}, {'P', -3}, {'S', -1}, {'T', -1}, {'W', -2}, {'Y', -2}, {'V', -1}}},
    {'Q', {{'A', -1}, {'R', 1}, {'N', 0}, {'D', 0}, {'C', -3}, {'Q', 5}, {'E', 2}, {'G', -2}, {'H', 0}, {'I', -3}, {'L', -2}, {'K', 1}, {'M', 0}, {'F', -3}, {'P', -1}, {'S', 0}, {'T', -1}, {'W', -2}, {'Y', -1}, {'V', -2}}},
    {'E', {{'A', -1}, {'R', 0}, {'N', 0}, {'D', 2}, {'C', -4}, {'Q', 2}, {'E', 5}, {'G', -2}, {'H', 0}, {'I', -3}, {'L', -3}, {'K', 1}, {'M', -2}, {'F', -3}, {'P', -1}, {'S', 0}, {'T', -1}, {'W', -3}, {'Y', -2}, {'V', -2}}},
    {'G', {{'A', 0}, {'R', -2}, {'N', 0}, {'D', -1}, {'C', -3}, {'Q', -2}, {'E', -2}, {'G', 6}, {'H', -2}, {'I', -4}, {'L', -4}, {'K', -2}, {'M', -3}, {'F', -3}, {'P', -2}, {'S', 0}, {'T', -2}, {'W', -2}, {'Y', -3}, {'V', -3}}},
    {'H', {{'A', -2}, {'R', 0}, {'N', 1}, {'D', -1}, {'C', -3}, {'Q', 0}, {'E', 0}, {'G', -2}, {'H', 8}, {'I', -3}, {'L', -3}, {'K', -1}, {'M', -2}, {'F', -1}, {'P', -2}, {'S', -1}, {'T', -2}, {'W', -2}, {'Y', 2}, {'V', -3}}},
    {'I', {{'A', -1}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', -1}, {'Q', -3}, {'E', -3}, {'G', -4}, {'H', -3}, {'I', 4}, {'L', 2}, {'K', -3}, {'M', 1}, {'F', 0}, {'P', -3}, {'S', -2}, {'T', -1}, {'W', -3}, {'Y', -1}, {'V', 3}}},
    {'L', {{'A', -1}, {'R', -2}, {'N', -3}, {'D', -4}, {'C', -1}, {'Q', -2}, {'E', -3}, {'G', -4}, {'H', -3}, {'I', 2}, {'L', 4}, {'K', -2}, {'M', 2}, {'F', 0}, {'P', -3}, {'S', -2}, {'T', -1}, {'W', -2}, {'Y', -1}, {'V', 1}}},
    {'K', {{'A', -1}, {'R', 2}, {'N', 0}, {'D', -1}, {'C', -3}, {'Q', 1}, {'E', 1}, {'G', -2}, {'H', -1}, {'I', -3}, {'L', -2}, {'K', 5}, {'M', -1}, {'F', -3}, {'P', -1}, {'S', 0}, {'T', -1}, {'W', -3}, {'Y', -2}, {'V', -2}}},
    {'M', {{'A', -1}, {'R', -1}, {'N', -2}, {'D', -3}, {'C', -1}, {'Q', 0}, {'E', -2}, {'G', -3}, {'H', -2}, {'I', 1}, {'L', 2}, {'K', -1}, {'M', 5}, {'F', 0}, {'P', -2}, {'S', -1}, {'T', -1}, {'W', -1}, {'Y', -1}, {'V', 1}}},
    {'F', {{'A', -2}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', -2}, {'Q', -3}, {'E', -3}, {'G', -3}, {'H', -1}, {'I', 0}, {'L', 0}, {'K', -3}, {'M', 0}, {'F', 6}, {'P', -4}, {'S', -2}, {'T', -2}, {'W', 1}, {'Y', 3}, {'V', -1}}},
    {'P', {{'A', -1}, {'R', -2}, {'N', -2}, {'D', -1}, {'C', -3}, {'Q', -1}, {'E', -1}, {'G', -2}, {'H', -2}, {'I', -3}, {'L', -3}, {'K', -1}, {'M', -2}, {'F', -4}, {'P', 7}, {'S', -1}, {'T', -1}, {'W', -4}, {'Y', -3}, {'V', -2}}},
    {'S', {{'A', 1}, {'R', -1}, {'N', 1}, {'D', 0}, {'C', -1}, {'Q', 0}, {'E', 0}, {'G', 0}, {'H', -1}, {'I', -2}, {'L', -2}, {'K', 0}, {'M', -1}, {'F', -2}, {'P', -1}, {'S', 4}, {'T', 1}, {'W', -3}, {'Y', -2}, {'V', -2}}},
    {'T', {{'A', 0}, {'R', -1}, {'N', 0}, {'D', -1}, {'C', -1}, {'Q', -1}, {'E', -1}, {'G', -2}, {'H', -2}, {'I', -1}, {'L', -1}, {'K', -1}, {'M', -1}, {'F', -2}, {'P', -1}, {'S', 1}, {'T', 5}, {'W', -2}, {'Y', -2}, {'V', 0}}},
    {'W', {{'A', -3}, {'R', -3}, {'N', -4}, {'D', -4}, {'C', -2}, {'Q', -2}, {'E', -3}, {'G', -2}, {'H', -2}, {'I', -3}, {'L', -2}, {'K', -3}, {'M', -1}, {'F', 1}, {'P', -4}, {'S', -3}, {'T', -2}, {'W', 11}, {'Y', 2}, {'V', -3}}},
    {'Y', {{'A', -2}, {'R', -2}, {'N', -2}, {'D', -3}, {'C', -2}, {'Q', -1}, {'E', -2}, {'G', -3}, {'H', 2}, {'I', -1}, {'L', -1}, {'K', -2}, {'M', -1}, {'F', 3}, {'P', -3}, {'S', -2}, {'T', -2}, {'W', 2}, {'Y', 7}, {'V', -1}}},
    {'V', {{'A', 0}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', -1}, {'Q', -2}, {'E', -2}, {'G', -3}, {'H', -3}, {'I', 3}, {'L', 1}, {'K', -2}, {'M', 1}, {'F', -1}, {'P', -2}, {'S', -2}, {'T', 0}, {'W', -3}, {'Y', -1}, {'V', 4}}}
};

void print_alignment(Alignment &alignment)
{
    int align_len = alignment.seq1.size();
    int align_len2 = alignment.seq2.size();

    if (align_len != align_len2){
        throw std::runtime_error("Alignment length mismatch");
    }

    int line_len = 80;
    cout << "Alignment score:" << alignment.score << endl;

    for (int line = 0; line < (align_len / line_len) + 1; line++){
        for (int i = line * line_len; (i < (line + 1) * line_len) && (i < align_len); i++){
            cout << alignment.seq1[i];
        }
        printf("\n");
        for (int i = line * line_len; (i < (line + 1) * line_len) && (i < align_len); i++){
            cout << ((alignment.seq1[i] == alignment.seq2[i]) ? '*' : ' ');
        }
        printf("\n");
        for (int i = line * line_len; (i < (line + 1) * line_len) && (i < align_len); i++){
            cout << alignment.seq2[i];
        }
        cout << endl;
    }
    return;
}

void check_seq(const std::string &seq){
    for (size_t i = 0; i < seq.size(); ++i){
        if ((seq[i] < 'A' || seq[i] > 'Z') && seq[i] != '-'){
            throw std::invalid_argument("Invalid character '" + std::string(1, seq[i]) + "' at position " + std::to_string(i) + " in sequence: " + seq);
        }
    }
}

std::unordered_map<char, std::unordered_map<char, int>> read_blosum(const std::string &file_path)
{
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + file_path);
    }

    std::unordered_map<char, std::unordered_map<char, int>> blosum_matrix;
    std::vector<char> amino_acids;
    std::string line;

    // Read the file line by line
    while (std::getline(file, line)) {
        // Skip comment lines or empty lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);

        // If this is the header line (amino acid labels), populate `amino_acids`
        if (line[0] == ' ' && amino_acids.empty()) {
            char aa;
            while (iss >> aa) {
                if (aa != 'X')
                    amino_acids.push_back(aa);
            }
            continue;
        }

        // Parse the matrix rows
        char row_label;
        iss >> row_label; // First column is the row label (amino acid)

        if (row_label == 'X') {
            continue; // Skip the 'X' row if present
        }

        std::unordered_map<char, int> row_values;
        for (size_t i = 0; i < amino_acids.size(); ++i) {
            int score;
            iss >> score;
            row_values[amino_acids[i]] = score;
        }

        blosum_matrix[row_label] = row_values;
    }

    return blosum_matrix;
}

int getScore(char a, char b, std::unordered_map<char, std::unordered_map<char, int>> blosum62)
{
    if (blosum62.find(a) != blosum62.end() && blosum62[a].find(b) != blosum62[a].end()) {
        return blosum62[a][b];
    }
    // cout << "Warning: Unknown characters '" << a << "' and '" << b << "' in BLOSUM matrix. Returning default penalty." << std::endl;
    return -4; // Default penalty for unknown characters
}

std::string removeGaps(const std::string &seq)
{
    std::string result = "";
    for (char c : seq) {
        if (c != '-') {
            result += c;
        }
    }
    return result;
}

Alignment sequence_align(const string seq1, const string seq2, const string matrix_file, int gapCost, int gapExtension)
{

    // Remove gaps from input sequences
    std::string seq1_nogap = removeGaps(seq1);
    std::string seq2_nogap = removeGaps(seq2);

    int len1 = seq1_nogap.length();
    int len2 = seq2_nogap.length();

    // Initialize the matrix
    std::vector<std::vector<int>> matrix(len1 + 1, std::vector<int>(len2 + 1, 0));
    std::vector<bool> prev_line(len2 + 1, false);

    std::unordered_map<char, std::unordered_map<char, int>> blosum62;
    if (matrix_file.empty()) {
        blosum62 = BLOSUM62;
        cout << "Using default BLOSUM62 matrix." << std::endl;
    } else {
        blosum62 = read_blosum(matrix_file);
    }

    // // Print blosum62 for debugging:
    // cout << "Loading BLOSUM62 matrix from: " << matrix_file << endl;
    // cout << "BLOSUM62 matrix loaded successfully." << endl;
    // // Print the BLOSUM62 matrix
    // cout << "BLOSUM62 matrix:" << endl;
    // for (const auto& row : blosum62) {
    //     cout << row.first << ": ";
    //     for (const auto& col : row.second) {
    //         cout << col.first << "=" << col.second << " ";
    //     }
    //     cout << endl;
    // }

    // Fill the matrix
    for (int i = 0; i < len1; i++) {
        bool prev = false; // insertion matrix[i, j - 1]
        for (int j = 0; j < len2; j++) {
            int choices[3];

            // Match/Mismatch
            choices[0] = matrix[i][j] + getScore(seq1_nogap[i], seq2_nogap[j], blosum62);

            // Delete (gap in seq2)
            choices[1] = matrix[i][j + 1] + (prev ? gapExtension : gapCost);

            // Insert (gap in seq1)
            choices[2] = matrix[i + 1][j] + (prev_line[j + 1] ? gapExtension : gapCost);

            // Find maximum
            int max_index = 0;
            for (int k = 1; k < 3; k++) {
                if (choices[k] >= choices[max_index]) {
                    max_index = k;
                }
            }

            matrix[i + 1][j + 1] = choices[max_index];
            prev_line[j + 1] = false;
            prev = false;

            if (max_index == 1) {
                prev = true;
            }
            else if (max_index == 2) {
                prev_line[j + 1] = true;
            }
        }
    }

    // cout << "Matrix filled successfully." << endl;
    // // Print the matrix for debugging
    // for (int i = 0; i <= len1; i++) {
    //     for (int j = 0; j <= len2; j++) {
    //         cout << matrix[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    // Find the maximum score (local alignment characteristic)
    int min_seq = std::min(len1, len2);
    int min_i = min_seq, min_j = min_seq;
    int max_score = matrix[min_seq][min_seq];

    for (int i = min_seq; i <= len1; i++) {
        for (int j = min_seq; j <= len2; j++) {
            if (matrix[i][j] > max_score) {
                max_score = matrix[i][j];
                min_i = i;
                min_j = j;
            }
        }
    }

    int i = min_i;
    int j = min_j;

    // Traceback and compute the alignment
    std::string align1 = "";
    std::string align2 = "";

    // Handle trailing sequences
    if (i != len1) {
        align2 = std::string(len1 - i, '-');
    }

    if (j != len2) {
        align1 = std::string(len2 - j, '-');
    }

    align1 += seq1_nogap.substr(i);
    align2 += seq2_nogap.substr(j);

    // Traceback through the matrix
    while (i != 0 && j != 0) {
        if (matrix[i][j] == matrix[i - 1][j - 1] + getScore(seq1_nogap[i - 1], seq2_nogap[j - 1], blosum62)) {
            align1 = seq1_nogap[i - 1] + align1;
            align2 = seq2_nogap[j - 1] + align2;
            i--;
            j--;
        } else if (matrix[i][j] == matrix[i - 1][j] + gapCost ||
                 matrix[i][j] == matrix[i - 1][j] + gapExtension) {
            align1 = seq1_nogap[i - 1] + align1;
            align2 = "-" + align2;
            i--;
        } else if (matrix[i][j] == matrix[i][j - 1] + gapCost ||
                 matrix[i][j] == matrix[i][j - 1] + gapExtension) {
            align1 = "-" + align1;
            align2 = seq2_nogap[j - 1] + align2;
            j--;
        } else {
            std::cerr << "Error in traceback" << std::endl;
            break;
        }
    }

    // Handle remaining prefixes
    align1 = seq1_nogap.substr(0, i) + align1;
    align2 = seq2_nogap.substr(0, j) + align2;

    if (i != 0) {
        align2 = std::string(i, '-') + align2;
    }
    else if (j != 0) {
        align1 = std::string(j, '-') + align1;
    }

    // Verify alignment lengths match
    assert(align1.length() == align2.length());

    Alignment alignment;
    alignment.seq1 = align1;
    alignment.seq2 = align2;
    alignment.score = max_score;

    return alignment;
}

void test(string seq_1, string seq_2)
{

    printf("Test\n");

    // Alignment alignment;

    Alignment alignment = sequence_align(seq_1, seq_2, "src/pdb_cpp/data/blosum62.txt", -11, -1);
    // printf ("Alignment:\n%s:\n%s:\n", alignment->seq1, alignment->seq2);

    print_alignment(alignment);
}

pair<vector<int>, vector<int>> get_common_atoms(
    const Coor &coor_1,
    const Coor &coor_2,
    const vector<string> &chain_1,
    const vector<string> &chain_2,
    const vector<string> &back_names,
    const string &matrix_file) {
    /**
     * Get atom selection in common for two Coor objects based on sequence alignment.
     *
     * Parameters:
     * - coor_1: First coordinate object
     * - coor_2: Second coordinate object
     * - chain_1: List of chains to consider in the first coordinate
     * - chain_2: List of chains to consider in the second coordinate
     * - back_names: List of backbone atom names
     *
     * Returns:
     * - pair of vectors containing indices for coor_1 and coor_2 respectively
     */

    // Helper function to join strings with spaces
    auto join_strings = [](const vector<string> &strings) -> string {
        if (strings.empty())
            return "";
        ostringstream oss;
        for (size_t i = 0; i < strings.size(); ++i) {
            if (i > 0)
                oss << " ";
            oss << strings[i];
        }
        return oss.str();
    };

    // Build selection strings
    string chain_1_str = join_strings(chain_1);
    string chain_2_str = join_strings(chain_2);
    string back_names_str = join_strings(back_names);

    string selection_1 = "chain " + chain_1_str +
                              " and protein and name " + back_names_str +
                              " and not altloc B C D E F";
    string selection_2 = "chain " + chain_2_str +
                              " and protein and name " + back_names_str +
                              " and not altloc B C D E F";

    // Select backbone atoms
    Coor coor_1_back = coor_1.select_atoms(selection_1);
    Coor coor_2_back = coor_2.select_atoms(selection_2);

    // Get sequences (assuming get_aa_sequences returns map<string, string>)
    auto sel_1_seq = coor_1_back.get_aa_sequences();
    auto sel_2_seq = coor_2_back.get_aa_sequences();

    // Get indices
    vector<int> sel_index_1 = coor_1.get_index_select(selection_1);
    vector<int> sel_index_2 = coor_2.get_index_select(selection_2);

    // Build concatenated sequences
    // Assuming get_aa_sequences() returns sequences in chain order matching the input chains
   string seq_1 = "";
    for (size_t i = 0; i < chain_1.size() && i < sel_1_seq.size(); ++i) {
        string chain_seq = sel_1_seq[i];
        // Remove gaps (replace "-" with "")
        chain_seq.erase(remove(chain_seq.begin(), chain_seq.end(), '-'), chain_seq.end());
        seq_1 += chain_seq;
    }
    
    string seq_2 = "";
    for (size_t i = 0; i < chain_2.size() && i < sel_2_seq.size(); ++i) {
        string chain_seq = sel_2_seq[i];
        // Remove gaps (replace "-" with "")
        chain_seq.erase(remove(chain_seq.begin(), chain_seq.end(), '-'), chain_seq.end());
        seq_2 += chain_seq;
    }

    // Assertions to check completeness
    assert(sel_index_1.size() == seq_1.length() * back_names.size() &&
           "Incomplete backbone atoms for first Coor object, you might consider using the remove_incomplete_residues method before.");
    assert(sel_index_2.size() == seq_2.length() * back_names.size() &&
           "Incomplete backbone atoms for second Coor object, you might consider using the remove_incomplete_residues method before.");

    // Perform sequence alignment
    Alignment alignment = sequence_align(seq_1, seq_2, matrix_file);
    string align_seq_1 = alignment.seq1;
    string align_seq_2 = alignment.seq2;

    // Build aligned selections
    vector<int> align_sel_1;
    vector<int> align_sel_2;

    size_t index_sel_1 = 0;
    size_t index_sel_2 = 0;
    size_t back_num = back_names.size();

    for (size_t i = 0; i < align_seq_1.length(); ++i) {
        if (align_seq_1[i] != '-' && align_seq_2[i] != '-') {
            // Both sequences have residues at this position
            for (size_t j = 0; j < back_num; ++j) {
                align_sel_1.push_back(sel_index_1[index_sel_1 + j]);
                align_sel_2.push_back(sel_index_2[index_sel_2 + j]);
            }
            index_sel_1 += back_num;
            index_sel_2 += back_num;
        } else if (align_seq_1[i] != '-') {
            // Only first sequence has residue
            index_sel_1 += back_num;
        } else {
            // Only second sequence has residue
            index_sel_2 += back_num;
        }
    }

    assert(align_sel_1.size() == align_sel_2.size() &&
           "Two selections don't have the same atom number");

    return make_pair(align_sel_1, align_sel_2);
}

/*
int main() {
    string seq_1 = "AQDMVSPPPPIADEPLTVNTGIYLIECYSLDDKAETFKVNAFLSLSWKDRRLAFDPV\
    RSGVRVKTYEPEAIWIPEIRFVNVENARDADVVDISVSPDGTVQYLERFSARVLSPLDFRRYPFDSQTLHIYLIVR\
    SVDTRNIVLAVDLEKVGKNDDVFLTGWDIESFTAVVKPANFALEDRLESKLDYQLRISRQYFSYIPNIILPMLFIL\
    FISWTAFWSTSYEANVTLVVSTLIAHIAFNILVETNLPKTPYMTYTGAIIFMIYLFYFVAVIEVTVQHYLKVESQP\
    ARAASITRASRIAFPVVFLLANIILAFLFFGF";
    string seq_2 = "MFALGIYLWETIVFFSLAASQQAAARKAASPMPPSEFLDKLMGKVSGYDARIRPNFK\
    GPPVNVTCNIFINSFGSIAETTMDYRVNIFLRQQWNDPRLAYSEYPDDSLDLDPSMLDSIWKPDLFFANEKGANFH\
    EVTTDNKLLRISKNGNVLYSIRITLVLACPMDLKNFPMDVQTCIMQLESFGYTMNDLIFEWDEKGAVQVADGLTLP\
    QFILKEEKDLRYCTKHYNTGKFTCIEARFHLERQMGYYLIQMYIPSLLIVILSWVSFWINMDAAPARVGLGITTVL\
    TMTTQSSGSRASLPKVSYVKAIDIWMAVCLLFVFSALLEYAAVNFIARQHKELLRFQRRRRHLKEDEAGDGRFSFA\
    AYGMGPACLQAKDGMAIKGNNNNAPTSTNPPEKTVEEMRKLFISRAKRIDTVSRVAFPLVFLIFNIFYWITYKIIR\
    SEDIHKQ";

    //string seq_1 = "AQDMVSPPXPIADEPLTVXSLSWKDRRL";
    //string seq_2 = "AQDMVSPPPPIADEPLTVNTGIYLIECYSLDDKAETFKVNAFLSLSWKDRRLAFDPV";



    cout << "START align" << endl;

    Alignment alignment;

    //alignment = sequence_align(seq_1, seq_2, -11, -1);


    //free(alignment->seq2);
    //free(alignment->seq1);

    int iter_num = 2;

    for (int i = 0; i < iter_num; i++) {
        cout << "Iter:" << i << endl;
        test(seq_1, seq_2);
    }


    return EXIT_SUCCESS;

}
*/
