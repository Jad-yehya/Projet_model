#include "tools.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define S 20
#define MIN_VAL (-1000)
#define MAX_VAL 1000

int main(int argc, char *argv[]) {

    srand(time(NULL));

    printf("========================================\n");
    printf("Starting the reconstruction error program \n");
    printf("========================================\n");

    FILE *f = fopen("float_recerr2.txt", "w");
    if (f == NULL) {
        fprintf(stderr, "Error opening file : float_recerr.txt\n");
        exit(1);
    }

    int nbtours = 2000;
    for(int k = 2; k < S; k++)
    {
        printf("k = %d \n", k);
        for (int i = 0; i < nbtours; ++i) {
            double **Q = (double **) malloc(k * sizeof(double *));
            double **R = (double **) malloc(k * sizeof(double *));
            double **A = (double **) malloc(k * sizeof(double *));
            for (int j = 0; j < k; j++) {
                A[j] = (double *) malloc(k * sizeof(double));
                Q[j] = (double *) malloc(k * sizeof(double));
                R[j] = (double *) malloc(k * sizeof(double));
            }

            generate_symmetric_matrix(&A, k, k, MIN_VAL, MAX_VAL);

            Givens3(A, k, k, &Q, &R);

            double **A_rec = (double **) malloc(k * sizeof(double *));
            for (int j = 0; j < k; j++) {
                A_rec[j] = (double *) malloc(k * sizeof(double));
            }

            multiplymatrices(Q, R, k, &A_rec);

            double rec_err = 0;
            for (int i = 0; i < k; ++i) {
                for (int j = 0; j < k; ++j) {
                    rec_err += fabs(A[i][j] - A_rec[i][j]);
                }
            }

            fprintf(f,"%d %f\n", k, rec_err);

            free_matrix(A, k);
            free_matrix(Q, k);
            free_matrix(R, k);
            free_matrix(A_rec, k);
 
        }   

    }   

    return 0;

}