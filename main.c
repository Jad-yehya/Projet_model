/**
* @authors CHAHINE Ezzedine, YEHYA Jad
*/
#include "tools.h"
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


    //Givens(A, 3, 2, Q, R);
    /*
    Givens2(A, 3, 2, &Q, &R);

    printf("In main : \n");


    printf("Q = \n");
    print_matrix(Q, 3, 3);

    printf("R = \n");
    print_matrix(R, 3, 2);
    */
    int nb_iter = 0;
    compute_eigenvalues(&A, 3, 3, &nb_iter);



    return 0;
}