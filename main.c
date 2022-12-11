/**
* @authors CHAHINE Ezzedine, YEHYA Jad
*/
#include "tools.h"
#include "tools_mpfr.h"


int main(int argc, char **argv) {
    double **A;
    double **Q;
    double **R;


    A = (double **) malloc(3 * sizeof(double *));
    for (int i = 0; i < 3; i++) {
        A[i] = (double *) malloc(2 * sizeof(double));
    }
    A[0][0] = 3;
    A[0][1] = -3;
    A[1][0] = 4;
    A[1][1] = -4;
    A[2][0] = 0;
    A[2][1] = 40;

    /*
    A[0][0] = 3;
    A[0][1] = 0;
    A[1][0] = 4;
    A[1][1] = 0;
    A[2][0] = 2.4;
    A[2][1] = 3.2;
    */
    printf("In main : A = \n");
    print_matrix(A, 3, 2);

    Q = (double **) malloc(3 * sizeof(double *));
    for (int i = 0; i < 3; i++) {
        Q[i] = (double *) malloc(3 * sizeof(double));
    }

    R = (double **) malloc(3 * sizeof(double *));
    for (int i = 0; i < 3; i++) {
        R[i] = (double *) malloc(2 * sizeof(double));
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
    compute_eigenvalues(&A,&Q, &R, 3, 2);

    return 0;
}