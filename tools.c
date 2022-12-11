//
// Created by Jad Yehya and Ezzedine Chahine on 24/10/2022.
//
// TODO Make a function read_file that reads a file and returns a matrix
// TODO Make a function write_file that writes a matrix to a file
#include "tools.h"

/**
 * @brief Create a matrix of size n*n
 * @param A
 * @param n
 */
void create_matrix(double **A, int n){
    A = (double **) malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        A[i] = (double *) malloc(n * sizeof(double));
    }
}

/**
 * @brief Swap two rows of a matrix
 * @param A : the matrix
 * @param i : the first row
 * @param j : the second row
 */
void swap_row(double **A, int i, int j) {
    double *temp = A[i];
    A[i] = A[j];
    A[j] = temp;
}


/**
 * Swap two columns of a matrix
 * @param A : the matrix
 * @param i : the first column
 * @param j : the second column
 * @param n : the number of rows
 */
void swap_columns(double **A, int i, int j, int n) {
    double temp;
    for (int k = 0; k < n; ++k) {
        temp = A[k][i];
        A[k][i] = A[k][j];
        A[k][j] = temp;
    }
}

/**
 * Print a matrix, we suppose that A is a square matrix
 * @param A : the matrix
 * @param m : number of rows
 * @param n : number of columns
 */
void print_matrix(double **A, int m, int n) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%f\t", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

/**
 * Free the matrix A
 * @param A : the matrix
 * @param n : the number of rows
 */
void free_matrix(double **A, int m) {
    for(int i = 0; i < m; i++) {
        free(A[i]);
    }
    free(A);
}

/**
 * Multiply two matrices
 * @param A : the first matrix
 * @param B : the second matrix
 * @param n : the size of the matrices
 * @return the result of the multiplication
 */
double **multiply_matrices(double **A, double **B, int n) {
    double **C = (double **) malloc(n * sizeof(double *));
    for (int i = 0; i < n; ++i) {
        C[i] = (double *) malloc(n * sizeof(double));
    }

    for (int i = 0; i < n; ++i) { // Rows of A
        for (int j = 0; j < n; ++j) { // Columns of B
            C[i][j] = 0;
            for (int k = 0; k < n; ++k) { // Columns of A
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (fabs(C[i][j]) < 1e-10) {
                C[i][j] = 0;
            }
        }
    }

    return C;
}

/**
 * Puts the indices of the maximum values of each row in an array in in and jn
 * @param A : Matrix
 * @param i : row index
 * @param j : column index
 * @param n : matrix size
 * @param in : (passage par pointeurs) max row index
 * @param jn : (passage par pointeurs) max column index
 */
void findGreatest(double **A, int i, int j, int n,  int *in, int *jn){
    double max = 0;
    for (int k = i; k < n-1; ++k) {
        for (int l = j; l < n-1; ++l) {
            if(fabs(A[k][l]) > max){
                max = fabs(A[k][l]);
                *in = k;
                *jn = l;
            }
        }
    }
}

/**
 * Generate the Givens matrix
 * @return
 */
void generate_G(double **G, int i, int j, int m, double** R){
    double c = R[j][j]/sqrt(R[j][j]*R[j][j] + R[i][j]*R[i][j]);
    double s = R[i][j]/sqrt(R[j][j]*R[j][j] + R[i][j]*R[i][j]);
    /*
    printf("----------------\n");
    printf("c = %f, s = %f\n", c, s);
    printf("----------------\n");
    */
    // Matrice G
    for (int k = 0; k < m; ++k) {
        if (k == i || k == j) {
            G[k][k] = c;
        } else {
            G[k][k] = 1;
        }
    }

    G[i][j] = -s;
    G[j][i] = s;



}

/**
 * Transpose the matrix M
 * @param M : Matrix
 * @param m : Number of rows
 * @param n : Number of columns
 * @param Mt : Transposed M
 */
void transpose_matrix(double ** M, int m, int n, double** Mt) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            Mt[j][i] = M[i][j];
        }
    }
}

/**
 * Copies A in B
 * @param A : Matrix to copy
 * @param B : Copy of A
 * @param m : Number of rows
 * @param n : Number of columns
 */
void copy(double **A, double **B, int m, int n){
    /*B = (double **) malloc(m * sizeof(double *));
    for (int i = 0; i < m; ++i) {
        B[i] = (double *) malloc(n * sizeof(double));
    }*/
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            B[i][j] = A[i][j];
        }
    }
}

/*void Givens(double **A, int m, int n, double** Q, double **R){
    // Matrice G
    double**G = (double **) malloc(m * sizeof(double *));
    for (int l = 0; l < m; ++l) {
        G[l] = (double *) malloc(m * sizeof(double));
    }

    // Matrice Gt
    double**Gt = (double **) malloc(m * sizeof(double *));
    for (int l = 0; l < m; ++l) {
        Gt[l] = (double *) malloc(m * sizeof(double));
    }

    // Initilazing Q
    for (int i = 0; i < m; ++i) {
        Q[i][i] = 1;
    }

    // Copying A in R
    copy(A, R, m, n);


    for (int j = 0; j < n; ++j) {
        for (int i = j+1; i < m; ++i) {
            // printf("Tour %d, %d\n", i, j);
            generate_G(G, i, j, m, R); //OK
            R = multiply_matrices(G, R, m);
            // print_matrix(G, m, m);
            transpose_matrix(G, m, m, Gt); //OK
            Q = multiply_matrices(Q, Gt, m);
            // Emptying G
            for (int k = 0; k < m; ++k) {
                for (int l = 0; l < m; ++l) {
                    G[k][l] = 0;
                }
            }

        }
    }
    print_matrix(Q, m, m);
    print_matrix(R, m, n);
}*/

void Givens2(double **A, int m, int n, double*** Q, double*** R){
    // Matrice G
    double**G = (double **) malloc(m * sizeof(double *));
    for (int l = 0; l < m; ++l) {
        G[l] = (double *) malloc(m * sizeof(double));
    }

    // Matrice Gt
    double**Gt = (double **) malloc(m * sizeof(double *));
    for (int l = 0; l < m; ++l) {
        Gt[l] = (double *) malloc(m * sizeof(double));
    }

    for (int i = 0; i < m; ++i) {
        (*Q)[i][i] = 1;
    }

    // Copying A in R
    copy(A, *R, m, n);

    for (int j = 0; j < n; ++j) {
        for (int i = j+1; i < m; ++i) {
            // printf("Tour %d, %d\n", i, j);
            generate_G(G, i, j, m, *R); //OK
            *R = multiply_matrices(G, *R, m);
            // print_matrix(G, m, m);
            transpose_matrix(G, m, m, Gt); //OK
            *Q = multiply_matrices(*Q, Gt, m);
            // Emptying G
            for (int k = 0; k < m; ++k) {
                for (int l = 0; l < m; ++l) {
                    G[k][l] = 0;
                }
            }

        }
    }
}

/**
 * Checks if the subdiagonal is 0
 * @param A : Matrix
 * @param m : Number of rows
 * @param n : Number of columns
 * @param epsilon : Precision
 * @return 0 if the subdiagonal is 0, 1 otherwise
 */
int thresh(double **A, int m, int n, double epsilon){
    int count = 0;
    for (int i = 1, j=0; i < m, j < n; ++i, ++j){
        if (fabs(A[i][j]) > epsilon){
            return 1;
        }
    }
    return 0;
}

/**
 * Recursive function that computes Q*A*Q^t
 * @param A Address of the matrix
 * @param m Number of rows
 * @param n Number of columns
 */
 // TODO : Remove Q and R from the parameters
 // TODO : Add a parameter epsilon
 // TODO : Add a parameter to count the number of iterations
 // TODO : Check errors of pointers
void compute_eigenvalues(double ***A, double ***Q, double ***R, int m, int n){
    printf("In compute_eigenvalues\n");
    if(thresh((*A), m, n, 1e-10) == 0){
        return;
    }
    /*
    double**Q = (double **) malloc(m * sizeof(double *));
    for (int l = 0; l < m; ++l) {
        Q[l] = (double *) malloc(m * sizeof(double));
    }

    double**R = (double **) malloc(m * sizeof(double *));
    for (int l = 0; l < m; ++l) {
        R[l] = (double *) malloc(n * sizeof(double));
    }
     */

    double **Qt = (double **) malloc(m * sizeof(double *));
    for (int l = 0; l < m; ++l) {
        Qt[l] = (double *) malloc(m * sizeof(double));
    }

    Givens2((*A), m, n, Q, R);
    transpose_matrix(*Q, m, m, Qt);
    printf("Q\n");
    print_matrix(*Q, m, m);
    printf("Q^t\n");
    print_matrix(Qt, m, m);

    (*A) = multiply_matrices(*Q, (*A), m);
    (*A) = multiply_matrices((*A), Qt, m);
    printf("A\n");
    print_matrix(*A, m, n);
    printf("\n");


    free_matrix(*Q, m);
    free_matrix(Qt, m);
    free_matrix(*R,  m);


    compute_eigenvalues(A, Q, R, m, n);
}
