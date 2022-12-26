//
// Created by Jad Yehya and Ezzedine Chahine on 24/10/2022.
//

#ifndef MODEL_PROJECT_TOOLS_MPFR_H
#define MODEL_PROJECT_TOOLS_MPFR_H

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
void swap_row_mpfr(mpfr_t **A, int row1, int row2, int n);
void swap_columns_mpfr(mpfr_t **A, int col1, int col2, int n);
void print_mpfr_matrix(mpfr_t **A, int n);
mpfr_t **multiply_matrices_mpfr(mpfr_t **A, mpfr_t **B, int n);
int check_equal_mpfr(mpfr_t **A, mpfr_t **B, int n);
#endif //MODEL_PROJECT_TOOLS_MPFR_H
