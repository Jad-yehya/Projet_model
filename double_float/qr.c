//
// Created by Jad Yehya on 24/12/2022.
//

#include "../tools.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define S 3
#define MIN_VAL (-10)
#define MAX_VAL 10

int main(int argc, char *argv[]) {

    srand(time(NULL));

    printf("========================================\n");
    printf("Starting the QR testing program \n");
    printf("========================================\n");

    FILE *f = fopen("QR_final.txt", "w");
    if (f == NULL) {
        fprintf(stderr, "Error opening file : QR.txt\n");
        exit(1);
    }

    int nbtours = 2000;
    clock_t start, end;
    double cpu_time_used;


    for (int k = 2; k < 10; k++) {
        printf("Size of the matrix : %d\n", k);
        for (int i = 0; i < 10; ++i) {
            double **Q = (double **) malloc(k * sizeof(double *));
            double **R = (double **) malloc(k * sizeof(double *));
            double **A = (double **) malloc(k * sizeof(double *));
            for (int j = 0; j < k; j++) {
                A[j] = (double *) malloc(k * sizeof(double));
                Q[j] = (double *) malloc(k * sizeof(double));
                R[j] = (double *) malloc(k * sizeof(double));
            }

            generate_symmetric_matrix(&A, k, k, MIN_VAL, MAX_VAL);

            start = clock();
            Givens3(A, k, k, &Q, &R);
            end = clock();

            cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
            fprintf(f,"%d %d %f\n", k, i, cpu_time_used);

            free_matrix(A, k);
            free_matrix(Q, k);
            free_matrix(R, k);

        }
    }

}