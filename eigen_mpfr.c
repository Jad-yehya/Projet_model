
#include "mpfr_impl.h"

#define M 4
#define N 4
#define MIN_VAL -10
#define MAX_VAL 10

int main(int argc, char const *argv[])
{
    printf("Starting the eigenvalues program\n");

    // Initialize the random number generator
    mpz_t seed;
    mpz_init(seed);
    gmp_randinit_default(RANDS);
    mpz_init_set_ui(seed, time(NULL));
    gmp_randseed(RANDS, seed);

    mpz_clear(seed);

    // Create a matrix of mpfr_t and initialize it
    mpfr_t **A = (mpfr_t **) malloc(M * sizeof(mpfr_t *));
    for (int i = 0; i < M; ++i) {
        A[i] = (mpfr_t *) malloc(N * sizeof(mpfr_t));
        for (int j = 0; j < N; ++j) {
            mpfr_init(A[i][j]);
        }
    }

    // Generate a symetric matrix
    generate_symetric_matrix_mpfr(&A, M, N, MIN_VAL, MAX_VAL);

    // Print the matrix
    printf("A = \n");
    print_mpfr_matrix(A, M, N);
    printf("\n\n\n");

    // Initialize threshold to 10^-10
    mpfr_t threshold;
    mpfr_init(threshold);
    mpfr_set_d(threshold, 1e-10, MPFR_RNDN);

    int nb_iter;
    compute_eigenvalues_mpfr(&A, M, N, &nb_iter, threshold);

    printf("Number of iterations: %d\n", nb_iter);

    // Frees
    free_mpfr_matrix(A, M, N);
    mpfr_clear(threshold);


    mpfr_free_cache();
    
    
    return 0;
}
