//
// Created by Jad Yehya and Ezzedine Chahine on 24/10/2022.
//

#include "tools_mpfr.h"


/**
 * Swaps two rows of a mpfr_t matrix
 * @param A : the matrix
 * @param row1 : the first row
 * @param row2 : the second row
 * @param n : size of the matrix
 */
void swap_row_mpfr(mpfr_t **A, int row1, int row2, int n) {
    mpfr_t temp;
    mpfr_init(temp);
    for (int i = 0; i < n; ++i) {
        mpfr_set(temp, A[row1][i], MPFR_RNDN);
        mpfr_set(A[row1][i], A[row2][i], MPFR_RNDN);
        mpfr_set(A[row2][i], temp, MPFR_RNDN);
    }
    mpfr_clear(temp);
}

/**
 * Swaps two columns of a mpfr_t matrix
 * @param A : the matrix
 * @param col1 : the first column
 * @param col2 : the second column
 * @param n : size of the matrix
 */
void swap_columns_mpfr(mpfr_t **A, int col1, int col2, int n) {
    mpfr_t temp;
    mpfr_init(temp);
    for (int i = 0; i < n; ++i) {
        mpfr_set(temp, A[i][col1], MPFR_RNDN);
        mpfr_set(A[i][col1], A[i][col2], MPFR_RNDN);
        mpfr_set(A[i][col2], temp, MPFR_RNDN);
    }
    mpfr_clear(temp);
}

/**
 * Print a matrix of mpfr_t
 * @param A : Matrix to print
 * @param n : Size of the matrix
 */
void print_mpfr_matrix(mpfr_t **A, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            mpfr_printf("%Rf\t", A[i][j]);
        }
        printf("\n");
    }
}

/**
 * Multiply two matrices of mpfr_t
 * @param A : the first matrix
 * @param B : the second matrix
 * @param n : the size of the matrices
 * @return the result of the multiplication
 */
mpfr_t **multiply_matrices_mpfr(mpfr_t **A, mpfr_t **B, int n) {
    mpfr_t **C = (mpfr_t **) malloc(n * sizeof(mpfr_t *));
    for (int i = 0; i < n; ++i) {
        C[i] = (mpfr_t *) malloc(n * sizeof(mpfr_t));
    }

    for (int i = 0; i < n; ++i) { // Rows of A
        for (int j = 0; j < n; ++j) { // Columns of B
            mpfr_init(C[i][j]);
            mpfr_set_d(C[i][j], 0, MPFR_RNDN);
            for (int k = 0; k < n; ++k) { // Columns of A
                mpfr_t temp;
                mpfr_init2(temp, 128);
                mpfr_mul(temp, A[i][k], B[k][j], MPFR_RNDN);
                mpfr_add(C[i][j], C[i][j], temp, MPFR_RNDN);
            }
        }
    }

    return C;
}

/**
 * Check if 2 mpfr_t matrices are equal
 * @param A : the first matrix
 * @param B : the second matrix
 * @param n : the size of the matrices
 * @return 0 if they are equal, 1 otherwise
 */
int check_equal_mpfr(mpfr_t **A, mpfr_t **B, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (mpfr_cmp(A[i][j], B[i][j]) != 0) {
                fprintf(stderr, "Matrices are not equal");
                return 1;
            }
        }
    }
    return 0;
}