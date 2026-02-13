#include "seq_align.h"
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>



int* seq_to_num (const char *seq)
{
    int seq_len = strlen(seq);
    // int *seq_num = malloc(seq_len * sizeof(int));
    int *seq_num = new int[seq_len];
    assert(seq_num != NULL);
    for (int i = 0; i < seq_len; i++)
    {
        switch (seq[i])
        {
            case 'A':
                seq_num[i] = 0;
                break;
            case 'R':
                seq_num[i] = 1;
                break;
            case 'N':
                seq_num[i] = 2;
                break;
            case 'D':
                seq_num[i] = 3;
                break;
            case 'C':
                seq_num[i] = 4;
                break;
            case 'Q':
                seq_num[i] = 5;
                break;
            case 'E':
                seq_num[i] = 6;
                break;
            case 'G':
                seq_num[i] = 7;
                break;
            case 'H':
                seq_num[i] = 8;
                break;
            case 'I':
                seq_num[i] = 9;
                break;
            case 'L':
                seq_num[i] = 10;
                break;
            case 'K':
                seq_num[i] = 11;
                break;
            case 'M':
                seq_num[i] = 12;
                break;
            case 'F':
                seq_num[i] = 13;
                break;
            case 'P':
                seq_num[i] = 14;
                break;
            case 'S':
                seq_num[i] = 15;
                break;
            case 'T':
                seq_num[i] = 16;
                break;
            case 'W':
                seq_num[i] = 17;
                break;
            case 'Y':
                seq_num[i] = 18;
                break;
            case 'V':
                seq_num[i] = 19;
                break;
            default:
                seq_num[i] = -1;
                break;
        }
    }

    return seq_num;
}

void print_alignment(Alignment *alignment)
{
    int align_len = strlen(alignment->seq1);
    int align_len2 = strlen(alignment->seq2);
    if (align_len != align_len2)
    {
        printf("Error: alignment length mismatch (%d, %d) !\n", align_len, align_len2);
        return;
    }
    int line_len = 80;
    printf("Alignment score: %d\n", alignment->score);

    for (int line = 0; line < (align_len / line_len) + 1; line++)
    {
        for (int i = line * line_len; (i < (line + 1) * line_len) && (i < align_len); i++)
        {
            printf("%c", alignment->seq1[i]);
        }
        printf("\n");
        for (int i = line * line_len; (i < (line + 1) * line_len) && (i < align_len); i++)
        {
            printf("%c",  (alignment->seq1[i] == alignment->seq2[i] ? '*' : ' '));
        }
        printf("\n");
        for (int i = line * line_len; (i < (line + 1) * line_len) && (i < align_len); i++)
        {
            printf("%c", alignment->seq2[i]);
        }
        printf("\n");

    }
    return ;
}

void read_matrix(const char *matrix_file, short int matrix[MATRIX_SIZE][MATRIX_SIZE])
{
    std::ifstream fp(matrix_file);
    if (!fp.is_open())
    {
        exit(EXIT_FAILURE);
    }

    std::string line;
    int i = 0;
    while (std::getline(fp, line) && i < MATRIX_SIZE)
    {
        if (line.empty() || line[0] == '#' || line[0] == ' ' || line[0] == '*')
            continue;

        std::istringstream iss(line);
        std::string row_label;
        iss >> row_label;

        for (int j = 0; j < MATRIX_SIZE; ++j)
        {
            int value = 0;
            if (!(iss >> value))
                break;
            matrix[i][j] = static_cast<short int>(value);
        }
        i++;
    }
}

// Alignment *align_test(const char *seq1, const char *seq2, const char *matrix_file, int GAP_COST, int GAP_EXT)
// {
//     // Alignment *alignment = malloc(sizeof(Alignment));
//     Alignment *alignment = new Alignment;
//     assert(alignment != NULL);
//     //alignment->seq1 = malloc((strlen(seq1) + 1) * sizeof(char));
//     alignment->seq1 = new char[strlen(seq1) + 1];
//     assert(alignment->seq1 != NULL);
//     //alignment->seq2 = malloc((strlen(seq2) + 1) * sizeof(char));
//     alignment->seq2 = new char[strlen(seq2) + 1];
//     assert(alignment->seq2 != NULL);
//     strcpy(alignment->seq1, seq1);
//     strcpy(alignment->seq2, seq2);
//     alignment->score = 0;
//     return alignment;
// }

void check_seq(const char *seq)
{
    int len = strlen(seq);
    for (int i = 0; i < len; i++)
    {
        if ((seq[i] < 'A' || seq[i] > 'Z') && seq[i] != '-')
        {
            printf("Invalid character '%c'|%d pos=%d\n", seq[i], seq[i], i);
            printf("Full sequence =  %s\n", seq);
            exit(EXIT_FAILURE);
        }
    }
}

//void align(const char *seq1, const char *seq2, const char *matrix_file, int GAP_COST, int GAP_EXT)


Alignment *seq_align(const char *seq1, const char *seq2, const char *matrix_file, int GAP_COST, int GAP_EXT)
{
    // std::cout << "Start seq_align "  << std::endl;
    short int subs_matrix[MATRIX_SIZE][MATRIX_SIZE];
    int seq1_len = strlen(seq1), seq2_len = strlen(seq2);
    // std::cout << "Start score matrix memory assignation"  << std::endl;
    // int score_matrix[seq1_len + 1][seq2_len + 1];
    // Allocate score matrix dynamically to avoid stack overflow with large sequences
    std::vector<std::vector<int>> score_matrix(seq1_len + 1, std::vector<int>(seq2_len + 1, 0));
    
    // std::cout << "End score matrix memory assignation"  << std::endl;
    std::vector<int> prev_score_line(seq2_len + 1, FALSE);
    int prev_score = FALSE;
    int i, j, counter;
    // Alignment *alignment = malloc(sizeof(Alignment));
    Alignment *alignment = new Alignment;



    assert(alignment != NULL);

    // clean input sequences
    check_seq(seq1);
    check_seq(seq2);

    // std::cout << "Reading substitution matrix from: " << matrix_file << std::endl;
    read_matrix(matrix_file, subs_matrix);

    // std::cout << "Substitution matrix loaded successfully." << std::endl;

    for (i = 0; i < seq1_len + 1; i++) score_matrix[i][0] = 0;
    for (i = 0; i < seq2_len + 1; i++) score_matrix[0][i] = 0;
    for (i = 0; i <= seq2_len; i++) prev_score_line[i] = FALSE;


    int *seq_1_num = seq_to_num(seq1);
    int *seq_2_num = seq_to_num(seq2);

    int match_score, insert_score, delete_score;

    // std::cout << "Starting matrix filling..." << std::endl;

    for (int i = 1; i <= seq1_len; i++)
    {
        prev_score = FALSE;
        for (int j = 1; j <= seq2_len; j++)
        {
            match_score = score_matrix[i - 1][j - 1] + subs_matrix[seq_1_num[i-1]][seq_2_num[j-1]];
            delete_score = score_matrix[i - 1][j] + (prev_score ? GAP_EXT : GAP_COST);
            insert_score = score_matrix[i][j - 1] + (prev_score_line[j] ? GAP_EXT : GAP_COST);

            if (match_score > insert_score && match_score > delete_score) {
                score_matrix[i][j] = match_score;
                prev_score_line[j] = FALSE;
                prev_score = FALSE;
            } else if (delete_score > insert_score) {
                score_matrix[i][j] = delete_score;
                prev_score_line[j] = FALSE;
                prev_score = TRUE;
            } else {
                score_matrix[i][j] = insert_score;
                prev_score_line[j] = TRUE;
                prev_score = FALSE;
            }
        }
    }

    // std::cout << "Matrix filled successfully." << std::endl;

    // Get max score index:
    int min_len = (seq1_len < seq2_len) ? seq1_len : seq2_len;
    int min_i=0, min_j=0, max_score = score_matrix[0][0];

    //int show_num = 15;
    //for (int i = 0; i < show_num; i++)
    //{
    //    for (int j = 0; j < show_num; j++)
    //    {
    //        printf ("%3d ", score_matrix[i][j]);
    //    }
    //    printf("\n");
    //}

    //int show_num = 10;
    //for (int i = seq1_len-show_num; i <= seq1_len; i++)
    //{
    //    for (int j = seq2_len-show_num; j <= seq2_len; j++)
    //    {const char *seq1, const char *seq2
    //        printf ("%3d ", score_matrix[i][j]);
    //    }
    //    printf("\n");
    //}
    //printf ("Get max \n");
    for (int i = min_len; i <= seq1_len; i++)
    {
        for (int j = min_len; j <= seq2_len; j++)
        {
            //printf ("%3d ", score_matrix[i][j]);
            if (score_matrix[i][j] > max_score) {
                max_score = score_matrix[i][j];
                min_i = i;
                min_j = j;
            }
        }
        //printf("\n");
    }

    //printf ("Max score: %d at %d, %d\n", max_score, min_i, min_j);

    alignment->score = max_score;

    // Traceback and compute the alignment using dynamic strings to avoid
    // fixed-buffer overflows.
    std::string rev_align_seq_1;
    std::string rev_align_seq_2;
    rev_align_seq_1.reserve(seq1_len + seq2_len + 2);
    rev_align_seq_2.reserve(seq1_len + seq2_len + 2);

    for (i = seq1_len - 1; i >= min_i; i--) {
        rev_align_seq_1.push_back(seq1[i]);
        rev_align_seq_2.push_back('-');
    }

    for (j = seq2_len - 1; j >= min_j; j--) {
        rev_align_seq_1.push_back('-');
        rev_align_seq_2.push_back(seq2[j]);
    }

    int trace_i = min_i - 1;
    int trace_j = min_j - 1;
    while (trace_i >= 0 && trace_j >= 0) {
        if (score_matrix[trace_i + 1][trace_j + 1] == score_matrix[trace_i][trace_j] + subs_matrix[seq_1_num[trace_i]][seq_2_num[trace_j]]) {
            rev_align_seq_1.push_back(seq1[trace_i--]);
            rev_align_seq_2.push_back(seq2[trace_j--]);
        }
        else if ((score_matrix[trace_i + 1][trace_j + 1] == score_matrix[trace_i][trace_j + 1] + GAP_COST) ||
                 (score_matrix[trace_i + 1][trace_j + 1] == score_matrix[trace_i][trace_j + 1] + GAP_EXT)) {
            rev_align_seq_1.push_back(seq1[trace_i--]);
            rev_align_seq_2.push_back('-');
        }
        else {
            rev_align_seq_1.push_back('-');
            rev_align_seq_2.push_back(seq2[trace_j--]);
        }
    }

    for (i = trace_i + 1; i > 0; i--) {
        rev_align_seq_1.push_back(seq1[i - 1]);
        rev_align_seq_2.push_back('-');
    }
    for (j = trace_j + 1; j > 0; j--) {
        rev_align_seq_1.push_back('-');
        rev_align_seq_2.push_back(seq2[j - 1]);
    }

    std::string final_align_seq_1(rev_align_seq_1.rbegin(), rev_align_seq_1.rend());
    std::string final_align_seq_2(rev_align_seq_2.rbegin(), rev_align_seq_2.rend());
    int align_len = static_cast<int>(final_align_seq_1.size());

    // Inverse the sequences
    //alignment->seq1 = calloc((align_len + 1), sizeof(char));
    alignment->seq1 = new char[align_len + 1];
    //alignment->seq2 = calloc((align_len + 1), sizeof(char));
    alignment->seq2 = new char[align_len + 1];
    assert(alignment->seq1 != NULL);
    assert(alignment->seq2 != NULL);
    //printf ("Alignment len align_len: %d\n", align_len + 1);

    for (int k = 0; k < align_len; ++k)
    {
        alignment->seq1[k] = final_align_seq_1[k];
        alignment->seq2[k] = final_align_seq_2[k];
    }
    alignment->seq1[align_len] = '\0';
    alignment->seq2[align_len] = '\0';

    //printf ("Alignment len %d\n", rev_i);

    delete[] seq_1_num;
    delete[] seq_2_num;
    check_seq(alignment->seq1);
    check_seq(alignment->seq2);

    //print_alignment(alignment);

    return alignment;
}

    // alignment = sequence_align(
    //     seq1=seq1,
    //     seq2=seq2,
    //     matrix_file='<package-resource blosum62.txt>',
    //     GAP_COST=gap_cost,
    //     GAP_EXT=gap_ext)

Alignment_cpp *sequence_align(
    const std::string &seq1, 
    const std::string &seq2, 
    const std::string &matrix_file, 
    int GAP_COST, int GAP_EXT) {
    // Convert C-style strings to C++ strings
    const char *c_seq1 = seq1.c_str();
    const char *c_seq2 = seq2.c_str();

    // Call the C function for sequence alignment
    Alignment *alignment = seq_align(c_seq1, c_seq2, matrix_file.c_str(), GAP_COST, GAP_EXT);

    // Create a C++ wrapper for the alignment result
    Alignment_cpp *alignment_cpp = new Alignment_cpp;
    alignment_cpp->seq1 = alignment->seq1;
    alignment_cpp->seq2 = alignment->seq2;
    alignment_cpp->score = alignment->score;

    // Free the C-style alignment structure
    free_align(alignment);

    return alignment_cpp;
}


void free_align(Alignment *align)
{
    if (align == NULL) {
        printf("Warning: alignment is NULL in free_align function !\n");
        return;
    }

    //printf("Freeing alignment seq1: %s\n", align->seq1);
    if (align->seq1 != NULL) {
        delete[] align->seq1;
        align->seq1 = NULL;
    } else {
        printf("Warning: seq1 is NULL\n");
    }
    //printf("Freeing alignment seq2: %s\n", align->seq2);
    if (align->seq2 != NULL) {
        delete[] align->seq2;
        align->seq2 = NULL;
    } else {
        printf("Warning: seq2 is NULL\n");
    }

    // Free memory for the alignment structure
    delete align;

}

// void test(char *seq_1, char *seq_2) {

//     printf("Test\n");

//     Alignment *alignment = NULL;
    
//     alignment = seq_align(seq_1, seq_2, "<matrix-file>", -11, -1);
//     printf ("Alignment:\n%s:\n%s:\n", alignment->seq1, alignment->seq2);

//     //print_alignment(alignment);

//     free_align(alignment);
// }

/*
int main(int argc, char *argv[])
{
    char *seq_1 = "AQDMVSPPPPIADEPLTVNTGIYLIECYSLDDKAETFKVNAFLSLSWKDRRLAFDPV\
RSGVRVKTYEPEAIWIPEIRFVNVENARDADVVDISVSPDGTVQYLERFSARVLSPLDFRRYPFDSQTLHIYLIVR\
SVDTRNIVLAVDLEKVGKNDDVFLTGWDIESFTAVVKPANFALEDRLESKLDYQLRISRQYFSYIPNIILPMLFIL\
FISWTAFWSTSYEANVTLVVSTLIAHIAFNILVETNLPKTPYMTYTGAIIFMIYLFYFVAVIEVTVQHYLKVESQP\
ARAASITRASRIAFPVVFLLANIILAFLFFGF";
    char *seq_2 = "MFALGIYLWETIVFFSLAASQQAAARKAASPMPPSEFLDKLMGKVSGYDARIRPNFK\
GPPVNVTCNIFINSFGSIAETTMDYRVNIFLRQQWNDPRLAYSEYPDDSLDLDPSMLDSIWKPDLFFANEKGANFH\
EVTTDNKLLRISKNGNVLYSIRITLVLACPMDLKNFPMDVQTCIMQLESFGYTMNDLIFEWDEKGAVQVADGLTLP\
QFILKEEKDLRYCTKHYNTGKFTCIEARFHLERQMGYYLIQMYIPSLLIVILSWVSFWINMDAAPARVGLGITTVL\
TMTTQSSGSRASLPKVSYVKAIDIWMAVCLLFVFSALLEYAAVNFIARQHKELLRFQRRRRHLKEDEAGDGRFSFA\
AYGMGPACLQAKDGMAIKGNNNNAPTSTNPPEKTVEEMRKLFISRAKRIDTVSRVAFPLVFLIFNIFYWITYKIIR\
SEDIHKQ";
    //char *seq_1 = "AQDMVSPPPPIADEPLTVNT";
    //char *seq_2 = "VSPPPPIADEP";

    short int sub_matrix[MATRIX_SIZE][MATRIX_SIZE];
    printf("Hello World!\n");
    
    printf ("READ MATRIX:\n");

    read_matrix("<matrix-file>", sub_matrix);

    printf ("MATRIX:\n");
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            printf("%2d ", sub_matrix[i][j]);
        }
        printf("\n");
    }

    //Alignment *alignment;

    //free(alignment->seq2);
    //free(alignment->seq1);

    int iter_num = 1000;

    for (int i = 0; i < iter_num; i++) {
        printf("Test %d\n", i);
        test(seq_1, seq_2);
    }


    return EXIT_SUCCESS;

}
*/