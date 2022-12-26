//
// Created by Jad Yehya and Ezzedine Chahine on 24/10/2022.
//
// TODO Make a function read_file that reads a file and returns a matrix
// TODO Make a function write_file that writes a matrix to a file
#include "tools.h"

/**
 * @brief Generate a random matrix of size m x n with values between a and b
 * (The matrix is allocated before calling this function)
 * @param A The matrix to print
 * @param m The number of rows
 * @param n The number of columns
 */
void generate_random_matrix(double ***A, int m, int n, int a, int b) {
    *A = (double **) malloc(m * sizeof(double *));
    for (int i = 0; i < m; i++) {
        (*A)[i] = (double *) malloc(n * sizeof(double));
    }
    double range = RAND_MAX / (b - a);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            (*A)[i][j] = (double) (rand() / range) + a;
        }
    }
}

/**
 * @brief Generate a symmetric matrix of size m x n with values between a and b
 * @param A The matrix to create
 * @param m The number of rows
 * @param n The number of columns
 * @param a The minimum value
 * @param b The maximum value
 */
void generate_symmetric_matrix(double ***A, int m, int n, int a, int b) {
    *A = (double **) malloc(m * sizeof(double *));
    for (int i = 0; i < m; i++) {
        (*A)[i] = (double *) malloc(n * sizeof(double));
    }
    double range = RAND_MAX / (b - a);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            (*A)[i][j] = ( ((double) rand()) / range) + a;
            (*A)[j][i] = (*A)[i][j];
        }
    }
}

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

void multiplymatrices(double **A, double **B, int n, double ***C) {
    for (int i = 0; i < n; ++i) { // Rows of A
        for (int j = 0; j < n; ++j) { // Columns of B
            (*C)[i][j] = 0;
            for (int k = 0; k < n; ++k) { // Columns of A
                (*C)[i][j] += A[i][k] * B[k][j];
            }
        }
    }
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
void transpose_matrix(double ** M, int m, int n, double*** Mt) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            (*Mt)[j][i] = M[i][j];
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
            transpose_matrix(G, m, m, &Gt); //OK
            *Q = multiply_matrices(*Q, Gt, m);
            // Emptying G
            for (int k = 0; k < m; ++k) {
                for (int l = 0; l < m; ++l) {
                    G[k][l] = 0;
                }
            }

        }
    }
    // FREES
    free_matrix(G, m);
    free_matrix(Gt, m);
}

void Givens3(double **A, int m, int n, double*** Q, double*** R){
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

    // R_tmp
    double**R_tmp = (double **) malloc(m * sizeof(double *));
    for (int l = 0; l < m; ++l) {
        R_tmp[l] = (double *) malloc(n * sizeof(double));
    }

    //Q_tmp
    double**Q_tmp = (double **) malloc(m * sizeof(double *));
    for (int l = 0; l < m; ++l) {
        Q_tmp[l] = (double *) malloc(m * sizeof(double));
    }


    // Copying A in R
    copy(A, *R, m, n);

    for (int j = 0; j < n; ++j) {
        for (int i = j+1; i < m; ++i) {
            generate_G(G, i, j, m, *R); //OK
            multiplymatrices(G, *R, m, &R_tmp);
            transpose_matrix(G, m, m, &Gt); //OK
            multiplymatrices(*Q, Gt, m, &Q_tmp);
            // Emptying G
            for (int k = 0; k < m; ++k) {
                for (int l = 0; l < m; ++l) {
                    G[k][l] = 0;
                }
            }
            // Copying R_tmp in R and Q_tmp in Q without using copy function
            for (int k = 0; k < m; ++k) {
                for (int l = 0; l < n; ++l) {
                    (*R)[k][l] = R_tmp[k][l];
                }
            }
            for (int k = 0; k < m; ++k) {
                for (int l = 0; l < m; ++l) {
                    (*Q)[k][l] = Q_tmp[k][l];
                }
            }
        }

    }
    // FREES
    free_matrix(G, m);
    free_matrix(Gt, m);
    free_matrix(R_tmp, m);
    free_matrix(Q_tmp, m);

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
    for (int i = 1; i < m; ++i) {
        if (fabs(A[i][i-1]) > epsilon) {
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
void compute_eigenvalues3(double ***A, int m, int n, long* nb_iter){
    if((thresh((*A), m, n, 1e-2) == 0) || (*nb_iter > 1000000)){
        printf("End of the recursion\n");
        printf("Number of iterations : %li\n", *nb_iter);
        printf("Eigenvalues : \n");
        for (int i = 0; i < m; ++i) {
            printf("%f\n", (*A)[i][i]);
        }
        return;
    }

    double**Q = (double **) malloc(m * sizeof(double *));
    for (int l = 0; l < m; ++l) {
        Q[l] = (double *) malloc(m * sizeof(double));
    }

    double**R = (double **) malloc(m * sizeof(double *));
    for (int l = 0; l < m; ++l) {
        R[l] = (double *) malloc(n * sizeof(double));
    }


    double **Qt = (double **) malloc(m * sizeof(double *));
    for (int l = 0; l < m; ++l) {
        Qt[l] = (double *) malloc(m * sizeof(double));
    }

    double **A2 = (double **) malloc(m * sizeof(double *));
    for (int l = 0; l < m; ++l) {
        A2[l] = (double *) malloc(n * sizeof(double));
    }

    double **A_temp = (double **) malloc(m * sizeof(double *));
    for (int l = 0; l < m; ++l) {
        A_temp[l] = (double *) malloc(n * sizeof(double));
    }

    // printf("In Eigenvalues3\n");

    Givens3((*A), m, n, &Q, &R);
    // printf("Givens done\n");
    transpose_matrix(Q, m, m, &Qt);
    // printf("Transpose done\n");

    //(*A) = multiply_matrices(Q, (*A), m);
    multiplymatrices(Q, (*A), m, &A2);
    //(*A) = multiply_matrices(A2, Qt, m);
    multiplymatrices(A2, Qt, m, &A_temp);

    // printf("Multiplication done\n");

    free_matrix(Q, m);
    free_matrix(Qt, m);
    free_matrix(R,  m);
    free_matrix(A2, m);

    // printf("Frees done\n");

    (*nb_iter)++;

    // A = A_temp without copy
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            (*A)[i][j] = A_temp[i][j];
        }
    }

    free_matrix(A_temp, m);


    compute_eigenvalues3(A, m, n, nb_iter);
}


/**
 * Iterative function that computes Q*A*Q^t
 * @param A
 * @param Q
 * @param R
 * @param m
 * @param n
 */
void compute_eigenvalues2(double **A, int m, int n)
{

    double**Q = (double **) malloc(m * sizeof(double *)); // Matrice Q taille m*m
    for (int l = 0; l < m; ++l)
    {
        Q[l] = (double *) malloc(m * sizeof(double));
    }

    double **Qt = (double **)malloc(m * sizeof(double *)); // Initialisation de Qt taille m*m
    for (int l = 0; l < m; ++l)
    {
        Qt[l] = (double *)malloc(m * sizeof(double));
    }

    double**R = (double **) malloc(m * sizeof(double *)); // Matrice R taille m*n
    for (int l = 0; l < m; ++l)
    {
        R[l] = (double *) malloc(n * sizeof(double));
    }

    double **A2 = (double **) malloc(m * sizeof(double *)); // Matrice A2 taille m*n
    for (int l = 0; l < m; ++l)
    {
        A2[l] = (double *) malloc(n * sizeof(double));
    }
    copy(A, A2, m, n);
    for (int i = 0; i < 100000; ++i)
    {
        // Threshold test
        if (thresh(A2, m, n, 1e-2) == 0)
        {
            printf("Threshold reached\n");
            printf("Number of iterations : %d\n", i);
            break;
        }
        Givens2(A2, m, n, &Q, &R);
        transpose_matrix(Q, m, m, &Qt);
        A2 = multiply_matrices(Q, A2, m);
        A2 = multiply_matrices(A2, Qt, m);
        //Emptying Q
        for (int k = 0; k < m; ++k)
        {
            for (int l = 0; l < m; ++l)
            {
                Q[k][l] = 0;
            }
        }
        //Emptying R
        for (int k = 0; k < m; ++k)
        {
            for (int l = 0; l < n; ++l)
            {
                R[k][l] = 0;
            }
        }
        //Emptying Qt
        for (int k = 0; k < m; ++k)
        {
            for (int l = 0; l < m; ++l)
            {
                Qt[k][l] = 0;
            }
        }
    }
    printf("A\n");
    print_matrix(A2, m, n);
    free_matrix(Q, m);
    free_matrix(Qt, m);
    free_matrix(R, m);
}