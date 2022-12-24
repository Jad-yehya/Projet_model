//
// Created by Jad Yehya on 23/12/2022.
//

#include "tools.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define S 4
#define MIN_VAL (-10)
#define MAX_VAL 10

int main(int argc, char **argv) {

    printf("========================================\n");
    printf("Starting the program with S = %d and MIN_VAL = %d and MAX_VAL = %d \n", S, MIN_VAL, MAX_VAL);
    printf("========================================\n");


    srand(time(NULL));

    // Generate a random matrix in for loop changing the size of the matrix and compute the eigenvalues
    /*for (int i = 2; i < 10; i++) {
        printf("Size of the matrix : %d\n", i);
        double** A = (double **) malloc(i * sizeof(double *));
        for (int j = 0; j < i; j++) {
            A[j] = (double *) malloc(i * sizeof(double));
        }
        generate_random_matrix(&A, i, i, 1, 10);
        printf("A = \n");
        print_matrix(A, i, i);


        int nb_iter = 0;

        compute_eigenvalues3(&A, i, i, &nb_iter);
        printf("Eigenvalues OK\n");
        free_matrix(A, i);
    }*/

    // Ouverture du fichier
    FILE *f = fopen("complexity_test_sym4.txt", "a");
    if (f == NULL) {
        fprintf(stderr, "Error opening file : complexity_test3.txt\n");
        exit(1);
    }

    long nb_iter;
    int nbtours = 500;
    clock_t start, end;
    double cpu_time_used;


    for (int i = 0; i < nbtours; i++) {
        if (i % 10 == 0) {
            printf("Tour %d\n", i);
        }
        nb_iter = 0;
        double** A = (double **) malloc(S * sizeof(double *));
        for (int j = 0; j < S; j++) {
            A[j] = (double *) malloc(S * sizeof(double));
        }
        // printf("Allocation OK\n");
        generate_symmetric_matrix(&A, S, S, MIN_VAL, MAX_VAL);

        // printf("Generation OK\n");
        // printf("A = \n");
        // print_matrix(A, S, S);

        start = clock();
        compute_eigenvalues3(&A, S, S, &nb_iter);
        end = clock();

        // printf("Eigenvalues OK\n");

        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        fprintf(f,"%li;%f\n",nb_iter, cpu_time_used);
        free_matrix(A, S);
    }

    fclose(f);






    return 0;
}
