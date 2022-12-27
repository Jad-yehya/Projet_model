#include "mpfr_impl.h"

#define M 20
#define N 4
#define MIN_VAL -1000
#define MAX_VAL 1000


int main() {
    // printf("Starting the reconstruction error estimation : \n");

    mpz_t seed;
    mpz_init(seed);
    gmp_randinit_default(RANDS);
    mpz_init_set_ui(seed, time(NULL));
    gmp_randseed(RANDS, seed);

    mpz_clear(seed);



    // Open a file to write the results
    FILE *f = fopen("reconstruction_error_quad.txt", "a");
    if (f == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }

    int nbtours = 10;
    for(int k = 2; k < M; k++) {
        // printf("k = %d \n", k);
        for (int i = 0; i < nbtours; ++i) {
            // Create a matrix of mpfr_t and initialize it
            mpfr_t **A = (mpfr_t **) malloc(k * sizeof(mpfr_t *));
            for (int i = 0; i < k; ++i) {
                A[i] = (mpfr_t *) malloc(k * sizeof(mpfr_t));
                for (int j = 0; j < k; ++j) {
                    mpfr_init(A[i][j]);
                }
            }

            // Generate a random matrix
            generate_symetric_matrix_mpfr(&A, k, k, MIN_VAL, MAX_VAL);

            // Create Q and R matrices
            mpfr_t **Q = (mpfr_t **) malloc(k * sizeof(mpfr_t *));
            mpfr_t **R = (mpfr_t **) malloc(k * sizeof(mpfr_t *));
            for (int i = 0; i < k; ++i) {
                Q[i] = (mpfr_t *) malloc(k * sizeof(mpfr_t));
                R[i] = (mpfr_t *) malloc(k * sizeof(mpfr_t));
                for (int j = 0; j < k; ++j) {
                    mpfr_init(Q[i][j]);
                    mpfr_init(R[i][j]);
                }
            }

            // Compute the QR decomposition
            Givens_mpfr(A, k, k, &Q, &R);

            // Compute the product of Q and R
            mpfr_t **QR = (mpfr_t **) malloc(k * sizeof(mpfr_t *));
            for (int i = 0; i < k; ++i) {
                QR[i] = (mpfr_t *) malloc(k * sizeof(mpfr_t));
                for (int j = 0; j < k; ++j) {
                    mpfr_init(QR[i][j]);
                }
            }

            // Compute the product of Q and R
            multiply_mpfr_matrices(Q, R, k, k, &QR);

            // Compute the reconstruction error
            mpfr_t error;
            mpfr_t tmp;
            mpfr_t quad;
            mpfr_init(error);
            mpfr_set_d(error, 0, MPFR_RNDN);
            mpfr_init(tmp);
            mpfr_init(quad);
            for (int i = 0; i < k; ++i) {
                mpfr_set_d(tmp, 0, MPFR_RNDN);
                for (int j = 0; j < k; ++j) {
                    mpfr_sub(tmp, A[i][j], QR[i][j], MPFR_RNDN);
                    mpfr_mul(quad, tmp, tmp, MPFR_RNDN);
                    mpfr_add(error, error, quad, MPFR_RNDN);

                }
            }

            // Write the error in the file
            fprintf(f, "%d ", k);
            mpfr_out_str(f, 10, 0, error, MPFR_RNDN);
            fprintf(f, "\n");

            // Free the memory
            free_mpfr_matrix(Q, k, k);
            free_mpfr_matrix(R, k, k);
            free_mpfr_matrix(QR, k, k);

            mpfr_clear(error);
            mpfr_clear(tmp);
            mpfr_clear(quad);
            free_mpfr_matrix(A, k, k);

            mpfr_free_cache();

        }
    }

    // Free the memory
    mpfr_free_cache();

    // Close the file
    fclose(f);



    return 0;
}


/**
 * - Problème : SEGFAULT MATRICE PLUS GRANDE QUE 5X5  OK
 * - R N'EST PAS UPPER TRIANGULAIRE
 * - IMPLEMENTER L ERREUR DE RECONSTRUCTION POUR LES FLOAT
 * - FAIRE LA COMPUTATION DES EIGENVALUES.
 */
