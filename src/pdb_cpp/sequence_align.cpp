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

#include "sequence_align.h"

using namespace std;

#define MATRIX_SIZE 23


vector<int> seq_to_num(const string &seq) {

    vector<int> seq_num(seq.size());

    for (size_t i = 0; i < seq.size(); i++) {
        switch (seq[i]) {
            case 'A': seq_num[i] = 0; break;
            case 'R': seq_num[i] = 1; break;
            case 'N': seq_num[i] = 2; break;
            case 'D': seq_num[i] = 3; break;
            case 'C': seq_num[i] = 4; break;
            case 'Q': seq_num[i] = 5; break;
            case 'E': seq_num[i] = 6; break;
            case 'G': seq_num[i] = 7; break;
            case 'H': seq_num[i] = 8; break;
            case 'I': seq_num[i] = 9; break;
            case 'L': seq_num[i] = 10; break;
            case 'K': seq_num[i] = 11; break;
            case 'M': seq_num[i] = 12; break;
            case 'F': seq_num[i] = 13; break;
            case 'P': seq_num[i] = 14; break;
            case 'S': seq_num[i] = 15; break;
            case 'T': seq_num[i] = 16; break;
            case 'W': seq_num[i] = 17; break;
            case 'Y': seq_num[i] = 18; break;
            case 'V': seq_num[i] = 19; break;
            case 'B': seq_num[i] = 20; break; // Asparagine or Aspartic acid
            case 'Z': seq_num[i] = 21; break; // Glutamine or Glutamic acid
            case 'X': seq_num[i] = 22; break; // Any amino acid
            default: seq_num[i] = 22; break;
        }
    }

    return seq_num;
}

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

void read_matrix(const string &matrix_file, vector<vector<short>> &matrix) {
    ifstream file(matrix_file);
    if (!file.is_open()) {
        throw runtime_error("Failed to open matrix file: " + matrix_file);
    }

    matrix.resize(MATRIX_SIZE, vector<short>(MATRIX_SIZE, 0));
    string line;
    int i = 0;

    while (getline(file, line) && i < MATRIX_SIZE) {
        if (line.empty() || line[0] == '#' || line[0] == ' ' || line[0] == '*') {
            continue; // Skip comments or empty lines
        }

        istringstream iss(line);
        string token;
        int j = 0;

        // Skip the first token (row label)
        iss >> token;

        // Read the matrix values
        while (iss >> token && j < MATRIX_SIZE) {
            matrix[i][j] = static_cast<short>(stoi(token));
            ++j;
        }
        ++i;
    }

    if (i < MATRIX_SIZE) {
        throw std::runtime_error("Matrix file does not contain enough rows.");
    }
}

void check_seq(const std::string &seq) {
    for (size_t i = 0; i < seq.size(); ++i) {
        if ((seq[i] < 'A' || seq[i] > 'Z') && seq[i] != '-') {
            throw std::invalid_argument("Invalid character '" + std::string(1, seq[i]) + "' at position " + std::to_string(i) + " in sequence: " + seq);
        }
    }
}

//void align(const char *seq1, const char *seq2, const char *matrix_file, int GAP_COST, int GAP_EXT)


Alignment sequence_align(const string seq1, const string seq2, const string matrix_file, int GAP_COST, int GAP_EXT) {
    vector<vector<short>> subs_matrix;
    int seq1_len = seq1.size(), seq2_len = seq2.size();
    int score_matrix [seq1_len + 1][seq2_len + 1];
    bool prev_score_line [seq2_len];
    bool prev_score;
    Alignment alignment;
    int i, j;

    // clean input sequences
    check_seq(seq1);
    check_seq(seq2);

    read_matrix(matrix_file, subs_matrix);

    for (i = 0; i < seq1_len + 1; i++) score_matrix[i][0] = 0;
    for (i = 0; i < seq2_len + 1; i++) score_matrix[0][i] = 0;
    for (i = 0; i < seq2_len; i++) prev_score_line[i] = false;


    vector<int> seq_1_num = seq_to_num(seq1);
    vector<int> seq_2_num = seq_to_num(seq2);

    int match, insert, del_flag;

    cout << "Start filling the score matrix:\n";
    cout << "Sequence 1: " << seq1 << endl;
    cout << "Sequence 2: " << seq2 << endl;
    cout << "Substitution matrix: " << matrix_file << endl;
    cout << "GAP_COST: " << GAP_COST << ", GAP_EXT: " << GAP_EXT << endl;
    cout << "Score matrix size: " << seq1_len + 1 << " x " << seq2_len + 1 << endl;
    cout << "Filling the score matrix:\n";
    // Fill the score matrix

    for (i = 1; i <= seq1_len; i++) {
        prev_score = false;
        //cout << i << ": " << endl;
        for (j = 1; j <= seq2_len; j++) {
            //cout << j << ": ";
            match = score_matrix[i - 1][j - 1] + subs_matrix[seq_1_num[i-1]][seq_2_num[j-1]];
            del_flag = score_matrix[i - 1][j] + (prev_score ? GAP_EXT : GAP_COST);
            insert = score_matrix[i][j - 1] + (prev_score_line[j] ? GAP_EXT : GAP_COST);

            if (match > insert && match > del_flag) {
                score_matrix[i][j] = match;
                prev_score_line[j] = false;
                prev_score = false;
            } else if (del_flag > insert) {
                score_matrix[i][j] = del_flag;
                prev_score_line[j] = false;
                prev_score = true;
            } else {
                score_matrix[i][j] = insert;
                prev_score_line[j] = true;
                prev_score = false;
            }
        }
        //cout << endl;
    }

    // cout << "Score matrix DONE:\n";
    // for (int i = 0; i <= seq1_len; i++) {
    //     for (int j = 0; j <= seq2_len; j++) {
    //         printf("%3d ", score_matrix[i][j]);
    //     }
    //     printf("\n");
    // }

    // Get max score index:
    int min_len = (seq1_len < seq2_len) ? seq1_len : seq2_len;
    int min_i=0, min_j=0, max_score = score_matrix[0][0];

    for (size_t i = min_len; i <= seq1_len; i++) {
        for (size_t j = min_len; j <= seq2_len; j++) {
            //printf ("%3d ", score_matrix[i][j]);
            if (score_matrix[i][j] > max_score) {
                max_score = score_matrix[i][j];
                min_i = i;
                min_j = j;
            }
        }
        //printf("\n");
    }

    printf ("Max score: %d at %d, %d\n", max_score, min_i, min_j);

    alignment.score = max_score;

    // Traceback and compute the alignment

    string align_seq_1="", align_seq_2="";

    for (i = seq1_len - 1; i >= min_i; i--) {
        align_seq_1 += seq1[i];
        align_seq_2 += '-';
    }

    for (j = seq2_len - 1; j >= min_j; j--) {
        align_seq_2 += seq2[j];
        align_seq_1 += '-';
    }


    printf ("Start matrix backtrack i=%d, j=%d, \n", i, j);
    cout << align_seq_1 << endl;
    cout << align_seq_2 << endl;

    do {
        //cout << "Backtrack i=" << i << ", j=" << j << endl;
        if (score_matrix[i+1][j+1] == score_matrix[i][j] + subs_matrix[seq_1_num[i]][seq_2_num[j]]) {
            align_seq_1 += seq1[i--];
            align_seq_2 += seq2[j--];
            //printf ("Match: %c, %c at i=%d, j=%d counter=%d\n", align_seq_1[counter], align_seq_2[counter], i, j, counter);
        }
        else if ((score_matrix[i+1][j+1] == score_matrix[i][j+1] + GAP_COST) ||
                 (score_matrix[i+1][j+1] == score_matrix[i][j+1] + GAP_EXT)) {
            align_seq_1 += seq1[i--];
            align_seq_2 += '-';
            //printf ("Insert: %c, %c at i=%d, j=%d counter=%d\n", align_seq_1[counter], align_seq_2[counter], i, j, counter);
        }
        else {
            align_seq_1 += '-';
            align_seq_2 += seq2[j--];
            //printf ("Delete: %c, %c at i=%d, j=%d counter=%d\n", align_seq_1[counter], align_seq_2[counter], i, j, counter);
        }
    } while (i >= 0 && j >= 0);


    //printf ("C Checking align sequences 1 len = %ld:\n", strlen(align_seq_1));
    //check_seq(align_seq_1);
    //printf ("C Checking align sequences 2 len = %ld:\n", strlen(align_seq_2));
    //check_seq(align_seq_2);

    for (i = i + 1; i > 0; i--) {
        align_seq_1 += seq1[i - 1];
        align_seq_2 += '-';
    }
    for (j = j + 1; j > 0; j--) {
        align_seq_2 += seq2[j - 1];
        align_seq_1 += '-';
    }

    //printf ("D Checking align sequences 1 len = %ld:\n", strlen(align_seq_1));
    //check_seq(align_seq_1);
    //printf ("D Checking align sequences 2 len = %ld:\n", strlen(align_seq_2));
    //check_seq(align_seq_2);
    // Get size of the alignment


    std::reverse(align_seq_1.begin(), align_seq_1.end());
    std::reverse(align_seq_2.begin(), align_seq_2.end());

    alignment.seq1 = align_seq_1;
    alignment.seq2 = align_seq_2;

    //print_alignment(alignment);

    return alignment;
}


void test(string seq_1, string seq_2) {

    printf("Test\n");

    //Alignment alignment;
    
    Alignment alignment = sequence_align(seq_1, seq_2, "src/pdb_cpp/data/blosum62.txt", -11, -1);
    //printf ("Alignment:\n%s:\n%s:\n", alignment->seq1, alignment->seq2);

    print_alignment(alignment);

}


int main() {
//     string seq_1 = "AQDMVSPPPPIADEPLTVNTGIYLIECYSLDDKAETFKVNAFLSLSWKDRRLAFDPV\
// RSGVRVKTYEPEAIWIPEIRFVNVENARDADVVDISVSPDGTVQYLERFSARVLSPLDFRRYPFDSQTLHIYLIVR\
// SVDTRNIVLAVDLEKVGKNDDVFLTGWDIESFTAVVKPANFALEDRLESKLDYQLRISRQYFSYIPNIILPMLFIL\
// FISWTAFWSTSYEANVTLVVSTLIAHIAFNILVETNLPKTPYMTYTGAIIFMIYLFYFVAVIEVTVQHYLKVESQP\
// ARAASITRASRIAFPVVFLLANIILAFLFFGF";
//     string seq_2 = "MFALGIYLWETIVFFSLAASQQAAARKAASPMPPSEFLDKLMGKVSGYDARIRPNFK\
// GPPVNVTCNIFINSFGSIAETTMDYRVNIFLRQQWNDPRLAYSEYPDDSLDLDPSMLDSIWKPDLFFANEKGANFH\
// EVTTDNKLLRISKNGNVLYSIRITLVLACPMDLKNFPMDVQTCIMQLESFGYTMNDLIFEWDEKGAVQVADGLTLP\
// QFILKEEKDLRYCTKHYNTGKFTCIEARFHLERQMGYYLIQMYIPSLLIVILSWVSFWINMDAAPARVGLGITTVL\
// TMTTQSSGSRASLPKVSYVKAIDIWMAVCLLFVFSALLEYAAVNFIARQHKELLRFQRRRRHLKEDEAGDGRFSFA\
// AYGMGPACLQAKDGMAIKGNNNNAPTSTNPPEKTVEEMRKLFISRAKRIDTVSRVAFPLVFLIFNIFYWITYKIIR\
// SEDIHKQ";

    string seq_1 = "AQDMVSPPXPIADEPLTVXSLSWKDRRL";
    string seq_2 = "AQDMVSPPPPIADEPLTVNTGIYLIECYSLDDKAETFKVNAFLSLSWKDRRLAFDPV";


    
    cout << "START align" << endl;

    //Alignment *alignment;

    //free(alignment->seq2);
    //free(alignment->seq1);

    int iter_num = 2;

    for (int i = 0; i < iter_num; i++) {
        cout << "Iter:" << i << endl;
        test(seq_1, seq_2);
    }


    return EXIT_SUCCESS;

}
