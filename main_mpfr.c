//
// Created by Jad Yehya on 25/12/2022.
//

#include "mpfr_impl.h"

#define M 4
#define N 4
#define MIN_VAL -10
#define MAX_VAL 10

int main()
{

    mpz_t seed;
    mpz_init(seed);
    gmp_randinit_default(RANDS);
    mpz_init_set_ui(seed, time(NULL));
    gmp_randseed(RANDS, seed);

    // Create a matrix of mpfr_t and initialize it
    mpfr_t **A = (mpfr_t **) malloc(M * sizeof(mpfr_t *));
    for (int i = 0; i < M; ++i) {
        A[i] = (mpfr_t *) malloc(N * sizeof(mpfr_t));
        for (int j = 0; j < N; ++j) {
            mpfr_init(A[i][j]);
        }
    }

    // Generate a random matrix
    generate_symetric_matrix_mpfr(&A, M, N, -10, 10);

    print_mpfr_matrix(A, M, N);
    free_mpfr_matrix(A, M, N);

    // Create a mpfr_t matrix of size 3x3 and intialize it with 1....9
    mpfr_t **B = (mpfr_t **) malloc(3 * sizeof(mpfr_t *));
    for (int i = 0; i < 3; ++i) {
        B[i] = (mpfr_t *) malloc(3 * sizeof(mpfr_t));
        for (int j = 0; j < 3; ++j) {
            mpfr_init_set_si(B[i][j], i * 3 + j + 1, MPFR_RNDN);
        }
    }

    printf("B = \n");
    print_mpfr_matrix(B, 3, 3);

    // Create a mpfr_t matrix of size 3x3 and intialize it with 3....11
    mpfr_t **C = (mpfr_t **) malloc(3 * sizeof(mpfr_t *));
    for (int i = 0; i < 3; ++i) {
        C[i] = (mpfr_t *) malloc(3 * sizeof(mpfr_t));
        for (int j = 0; j < 3; ++j) {
            mpfr_init_set_si(C[i][j], i * 3 + j + 3, MPFR_RNDN);
        }
    }

    printf("C = \n");
    print_mpfr_matrix(C, 3, 3);

    // Create a mpfr_t matrix of size 3x3 and intialize it with 0
    mpfr_t **D = (mpfr_t **) malloc(3 * sizeof(mpfr_t *));
    for (int i = 0; i < 3; ++i) {
        D[i] = (mpfr_t *) malloc(3 * sizeof(mpfr_t));
        for (int j = 0; j < 3; ++j) {
            mpfr_init(D[i][j]);
        }
    }

    // Copy B into D
    copy_mpfr_matrix(B, 3, 3, &D);
    

    // Free the memory
    free_mpfr_matrix(B, 3, 3);
    free_mpfr_matrix(C, 3, 3);
    free_mpfr_matrix(D, 3, 3);


    return 0;
}