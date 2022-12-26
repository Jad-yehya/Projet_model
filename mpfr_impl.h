
#ifndef MODEL_PROJET_MPFR_IMPL_H_
# define MODEL_PROJET_MPFR_IMPL_H_

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include <gmp.h>
# include <mpfr.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

// Seed for the random number generator.
// That way we only seed once.
gmp_randstate_t RANDS;



void generate_random_matrix_mpfr(mpfr_t *** matrix, int m, int n, int a, int b);
void generate_symetric_matrix_mpfr(mpfr_t *** matrix, int m, int n, int a, int b);
void print_mpfr_matrix(mpfr_t ** matrix, int m, int n);
void multiply_mpfr_matrices(mpfr_t **A, mpfr_t **B, int m, int n, mpfr_t ***C);
void generate_givens_matrix(mpfr_t *** G, int m, int i, int j, mpfr_t **R);
void free_mpfr_matrix(mpfr_t ** matrix, int m, int n);
void transpose_mpfr_matrix(mpfr_t **M, int m, int n, mpfr_t ***Mt);
void copy_mpfr_matrix(mpfr_t **A, int m, int n, mpfr_t ***B);
void Givens_mpfr(mpfr_t ** A, int m, int n, mpfr_t*** Q, mpfr_t*** R);
int thresh_mpfr(mpfr_t ** A, int m, int n, mpfr_t threshold);
void compute_eigenvalues_mpfr(mpfr_t*** A, int m, int n, int* nb_iter, mpfr_t threshold);

#endif //MODEL_PROJET_MPFR_IMPL_H_