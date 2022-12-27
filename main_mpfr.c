//
// Created by Jad Yehya on 25/12/2022.
//

#include "mpfr_impl.h"

#define M 4
#define N 4
#define MIN_VAL -1000
#define MAX_VAL 1000

int main()
{
    printf("Hello, World!\n");

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
    generate_symetric_matrix_mpfr(&A, M, N, MIN_VAL, MAX_VAL);

    print_mpfr_matrix(A, M, N);

    printf("MATRICE OK\n");

    // Create Q and R matrices
    mpfr_t **Q = (mpfr_t **) malloc(M * sizeof(mpfr_t *));
    mpfr_t **R = (mpfr_t **) malloc(M * sizeof(mpfr_t *));
    for (int i = 0; i < M; ++i) {
        Q[i] = (mpfr_t *) malloc(N * sizeof(mpfr_t));
        R[i] = (mpfr_t *) malloc(N * sizeof(mpfr_t));
        for (int j = 0; j < N; ++j) {
            mpfr_init(Q[i][j]);
            mpfr_init(R[i][j]);
        }
    }

    printf("Entering Givens_mpfr \n");
    // Compute the QR decomposition
    Givens_mpfr(A, M, N, &Q, &R);

    // Print the results
    printf("Q = \n");
    print_mpfr_matrix(Q, M, N);

    printf("R = \n");
    print_mpfr_matrix(R, M, N);

    // Compute the product of Q and R
    mpfr_t **QR = (mpfr_t **) malloc(M * sizeof(mpfr_t *));
    for (int i = 0; i < M; ++i) {
        QR[i] = (mpfr_t *) malloc(N * sizeof(mpfr_t));
        for (int j = 0; j < N; ++j) {
            mpfr_init(QR[i][j]);
        }
    }
    multiply_mpfr_matrices(Q, R, M, N, &QR);

    // Print the result
    printf("QR = \n");
    print_mpfr_matrix(QR, M, N);


    // Free the memory
    free_mpfr_matrix(A, M, N);
    free_mpfr_matrix(Q, M, N);
    free_mpfr_matrix(R, M, N);
    free_mpfr_matrix(QR, M, N);


    return 0;
}