/**
* @authors CHAHINE Ezzedine, YEHYA Jad
*/
#include "../tools.h"
//#include "tools_mpfr.h"


int main(int argc, char **argv) {
    double **A;
    double **Q;
    double **R;

    int m = 3;
    int n = 3;


    A = (double **) malloc(m * sizeof(double *));
    for (int i = 0; i < 3; i++) {
        A[i] = (double *) malloc(n * sizeof(double));
    }
    A[0][0] = 4;
    A[0][1] = 3;
    A[0][2] = 9;
    A[1][0] = 2;
    A[1][1] = 8;
    A[1][2] = 5;
    A[2][0] = 3;
    A[2][1] = 1;
    A[2][2] = 2;

    printf("In main : A = \n");
    print_matrix(A, m, m);

    Q = (double **) malloc(m * sizeof(double *));
    for (int i = 0; i < 3; i++) {
        Q[i] = (double *) malloc(m * sizeof(double));
    }

    R = (double **) malloc(m * sizeof(double *));
    for (int i = 0; i < 3; i++) {
        R[i] = (double *) malloc(m * sizeof(double));
    }

    long nb_iter = 0;
    // compute_eigenvalues3(&A, 3, 3, &nb_iter);
    Givens3(A, m, n, &Q, &R);

    printf("In main : Q = \n");
    print_matrix(Q, m, m);
    printf("In main : R = \n");
    print_matrix(R, m, m);


    free_matrix(A, m);
    free_matrix(Q, m);
    free_matrix(R, m);


    return 0;
}