#ifndef _ALIGN_H
#define _ALIGN_H

// Include necessary headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// Constants
#define MATRIX_SIZE 20
#define TRUE 1
#define FALSE 0

// Data structures
typedef struct {
    char *seq1;
    char *seq2;
    int score;
} Alignment;

// Function declarations

/**
 * Convert amino acid sequence to numerical representation
 * @param seq Input amino acid sequence as string
 * @return Pointer to array of integers representing amino acids
 */
int* seq_to_num(const char *seq);

/**
 * Print alignment in formatted output
 * @param alignment Pointer to Alignment structure to print
 */
void print_alignment(Alignment *alignment);

/**
 * Read substitution matrix from file
 * @param matrix_file Path to matrix file
 * @param matrix 2D array to store the substitution matrix
 */
void read_matrix(const char *matrix_file, short int matrix[MATRIX_SIZE][MATRIX_SIZE]);

/**
 * Create test alignment (simple copy of sequences)
 * @param seq1 First sequence
 * @param seq2 Second sequence
 * @param matrix_file Path to substitution matrix file
 * @param GAP_COST Gap opening cost
 * @param GAP_EXT Gap extension cost
 * @return Pointer to new Alignment structure
 */
Alignment *align_test(const char *seq1, const char *seq2, const char *matrix_file, int GAP_COST, int GAP_EXT);

/**
 * Perform sequence alignment using dynamic programming
 * @param seq1 First sequence to align
 * @param seq2 Second sequence to align
 * @param matrix_file Path to substitution matrix file
 * @param GAP_COST Gap opening cost
 * @param GAP_EXT Gap extension cost
 * @return Pointer to Alignment structure containing aligned sequences
 */
Alignment *seq_align(const char *seq1, const char *seq2, const char *matrix_file, int GAP_COST=-11, int GAP_EXT=-1);

/**
 * Validate sequence contains only valid amino acid characters
 * @param seq Sequence to check
 */
void check_seq(const char *seq);

/**
 * Free memory allocated for Alignment structure
 * @param align Pointer to Alignment structure to free
 */
void free_align(Alignment *align);

/**
 * Test function for alignment
 * @param seq_1 First sequence
 * @param seq_2 Second sequence
 */
void test(char *seq_1, char *seq_2);

#endif // _ALIGN_H
