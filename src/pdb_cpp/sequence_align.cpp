#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <memory>
#include <algorithm>
#include <cassert>
#include <unordered_map>

#include "sequence_align.h"

using namespace std;

#define MATRIX_SIZE 23


void print_alignment(Alignment &alignment) {
    int align_len = alignment.seq1.size();
    int align_len2 = alignment.seq2.size();

    if (align_len != align_len2) {
        throw std::runtime_error("Alignment length mismatch");
    }

    int line_len = 80;
    cout << "Alignment score:" << alignment.score << endl;

    for (int line = 0; line < (align_len / line_len) + 1; line++) {
        for (int i = line * line_len; (i < (line + 1) * line_len) && (i < align_len); i++) {
           cout << alignment.seq1[i];
        }
        printf("\n");
        for (int i = line * line_len; (i < (line + 1) * line_len) && (i < align_len); i++) {
            cout << ((alignment.seq1[i] == alignment.seq2[i]) ? '*' : ' ');
        }
        printf("\n");
        for (int i = line * line_len; (i < (line + 1) * line_len) && (i < align_len); i++) {
            cout << alignment.seq2[i];
        }
        cout << endl;
    }
    return ;
}

void check_seq(const std::string &seq) {
    for (size_t i = 0; i < seq.size(); ++i) {
        if ((seq[i] < 'A' || seq[i] > 'Z') && seq[i] != '-') {
            throw std::invalid_argument("Invalid character '" + std::string(1, seq[i]) + "' at position " + std::to_string(i) + " in sequence: " + seq);
        }
    }
}

std::unordered_map<char, std::unordered_map<char, int>> read_blosum(const std::string &file_path) {
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
                if (aa != 'X') amino_acids.push_back(aa);
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


int getScore(char a, char b, std::unordered_map<char, std::unordered_map<char, int>> blosum62) {
    if (blosum62.find(a) != blosum62.end() && blosum62[a].find(b) != blosum62[a].end()) {
        return blosum62[a][b];
    }
    //cout << "Warning: Unknown characters '" << a << "' and '" << b << "' in BLOSUM matrix. Returning default penalty." << std::endl;
    return -4; // Default penalty for unknown characters
}


std::string removeGaps(const std::string& seq) {
    std::string result = "";
    for (char c : seq) {
        if (c != '-') {
            result += c;
        }
    }
    return result;
}


Alignment sequence_align(const string seq1, const string seq2, const string matrix_file, int gapCost, int gapExtension) {

    // Remove gaps from input sequences
    std::string seq1_nogap = removeGaps(seq1);
    std::string seq2_nogap = removeGaps(seq2);

    int len1 = seq1_nogap.length();
    int len2 = seq2_nogap.length();

    // Initialize the matrix
    std::vector<std::vector<int>> matrix(len1 + 1, std::vector<int>(len2 + 1, 0));
    std::vector<bool> prev_line(len2 + 1, false);

    std::unordered_map<char, std::unordered_map<char, int>> blosum62 = read_blosum(matrix_file);
    
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
            } else if (max_index == 2) {
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
    } else if (j != 0) {
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

void test(string seq_1, string seq_2) {

    printf("Test\n");

    //Alignment alignment;
    
    Alignment alignment = sequence_align(seq_1, seq_2, "src/pdb_cpp/data/blosum62.txt", -11, -1);
    //printf ("Alignment:\n%s:\n%s:\n", alignment->seq1, alignment->seq2);

    print_alignment(alignment);

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

    string seq_1 = "AQDMVSPPXPIADEPLTVXSLSWKDRRL";
    string seq_2 = "AQDMVSPPPPIADEPLTVNTGIYLIECYSLDDKAETFKVNAFLSLSWKDRRLAFDPV";


    
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