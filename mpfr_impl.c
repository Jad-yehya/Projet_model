/**
 * @file mpfr_impl.c
 * @author YEHYA Jad
 * @brief Implementation of the project using the MPFR library
 * @version 0.1
 * @date 2022-12-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "mpfr_impl.h"

/***
 * PROBLEME SEGFAULT AVANT DE RENTRER DANS GENERATE GIVENS ROTATION MATRIX
 * TODO :
 * - RESOUDRE LE PROBLEME DE SEGFAULT
 * - VERIFIER QUE LES CALCULS SONT BONS
 * - CONTINUER A IMPLEMENTER LE RESTE 
 * - !!!!!!!!!!!!!!!!!!!! VERIFIER QUE G EST BONNE POUR LA PARTIE FLOATING POINT !!!!!!!!!
 * 
*/

/**
 * @brief Generate a random matrix of size m x n and store it in the matrix parameter.
 * The matrix is allocated before calling this function.
 * @param matrix 
 * @param m 
 * @param n 
 * @param a 
 * @param b 
 */
void generate_random_matrix_mpfr(mpfr_t *** matrix, int m, int n, int a, int b)
{
    mpfr_t r, min, max, range;
    mpfr_init2(r, 128);
    mpfr_init2(min, 128);
    mpfr_init2(max, 128);
    mpfr_init2(range, 128);

    mpfr_set_d(min, a, MPFR_RNDN);
    mpfr_set_d(max, b, MPFR_RNDN);
    mpfr_sub(range, max, min, MPFR_RNDN);

    for (int i = 0; i < m; i++)
    {
        for (int j = i+1; j < n; j++)
        {
            mpfr_urandom(r, RANDS, MPFR_RNDN);
            mpfr_mul(r, r, range, MPFR_RNDN);
            mpfr_add(r, r, min, MPFR_RNDN);

            mpfr_set((*matrix)[i][j], r, MPFR_RNDN);
        }
    }
    mpfr_clear(r);
    mpfr_clear(min);
    mpfr_clear(max);
    mpfr_clear(range);
}


/**
 * @brief Generate a random symetric matrix of size m x n and store it in the matrix parameter.
 * The matrix is allocated before calling this function.
 * @param matrix 
 * @param m 
 * @param n 
 * @param a 
 * @param b 
 */
void generate_symetric_matrix_mpfr(mpfr_t *** matrix, int m, int n, int a, int b)
{
    mpfr_t r, min, max, range;
    mpfr_init2(r, 128);
    mpfr_init2(min, 128);
    mpfr_init2(max, 128);
    mpfr_init2(range, 128);

    mpfr_set_d(min, a, MPFR_RNDN);
    mpfr_set_d(max, b, MPFR_RNDN);
    mpfr_sub(range, max, min, MPFR_RNDN);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            mpfr_urandom(r, RANDS, MPFR_RNDN);
            mpfr_mul(r, r, range, MPFR_RNDN);
            mpfr_add(r, r, min, MPFR_RNDN);

            mpfr_set((*matrix)[i][j], r, MPFR_RNDN);
            mpfr_set((*matrix)[j][i], r, MPFR_RNDN);
        }
    }
    mpfr_clear(r);
    mpfr_clear(min);
    mpfr_clear(max);
    mpfr_clear(range);

}


/**
 * @brief Print a matrix of mpfr_t
 * @param matrix 
 * @param m 
 */
void print_mpfr_matrix(mpfr_t ** matrix, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            mpfr_printf("%.16R*f\t", MPFR_RNDN, matrix[i][j]);
        }
        printf("\n");
    }
}


/**
 * @brief Multiply two matrices A and B and store the result in C.  
 * 
 * @param A First matrix
 * @param B Second matrix
 * @param m Number of rows of A
 * @param n Number of columns of A and rows of B
 * @param C Result matrix
 */
void multiply_mpfr_matrices(mpfr_t **A, mpfr_t **B, int m, int n, mpfr_t ***C)
{
    mpfr_t sum, tmp;
    mpfr_init2(sum, 128);
    mpfr_init2(tmp, 128);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            mpfr_set_d(sum, 0, MPFR_RNDN);
            for (int k = 0; k < n; k++)
            {
                mpfr_mul(tmp, A[i][k], B[k][j], MPFR_RNDN);
                mpfr_add(sum, sum, tmp, MPFR_RNDN);
            }
            mpfr_set((*C)[i][j], sum, MPFR_RNDN);
        }
    }
    mpfr_clear(sum);
    mpfr_clear(tmp);
}

/**
 * @brief Generate a Givens rotation matrix G(i,j). The matrix is allocated before calling this function.
 * 
 * @param G Result matrix
 * @param m Size of the matrix (m x m)
 * @param i Givens index
 * @param j Givens index
 * @param R Matrix to be decomposed
 */
void generate_givens_matrix(mpfr_t *** G, int m, int i, int j, mpfr_t **R)
{
    printf("DANS GENERATE GIVENS ROTATION MATRIX");

    // Making G an identity matrix
    for(int i=0; i < m; i++)
    {
        for(int j=0; j < m; j++)
        {
            if(i == j)
            {
                mpfr_set_d((*G)[i][j], 1, MPFR_RNDN);
            }
            else
            {
                mpfr_set_d((*G)[i][j], 0, MPFR_RNDN);
            }
        }
    }

    printf("IDENTITY MATRIX OK");

    mpfr_t c, s, tmp, tmp2;
    mpfr_init2(c, 128);
    mpfr_init2(s, 128);
    mpfr_init2(tmp, 128);

    mpfr_mul(tmp, R[j][j], R[j][j], MPFR_RNDN);
    mpfr_mul(tmp2, R[i][j], R[i][j], MPFR_RNDN);
    mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);

    mpfr_div(c, R[j][j], tmp, MPFR_RNDN);
    mpfr_div(s, R[i][j], tmp, MPFR_RNDN);

    printf("C AND S OK");

    mpfr_set((*G)[j][i], s, MPFR_RNDN);
    mpfr_neg(s, s, MPFR_RNDN);
    mpfr_set((*G)[i][j], s, MPFR_RNDN);

    for (int k = 0; k < m; k++)
    {
        if (k == i || k == j)
        {
            mpfr_set((*G)[k][k], c, MPFR_RNDN);
        }
    }
    mpfr_clear(c);
    mpfr_clear(s);
    mpfr_clear(tmp);
    mpfr_clear(tmp2);
}


/**
 * @brief Free the memory allocated for a matrix of mpfr_t
 * 
 * @param matrix Matrix to free
 * @param m Number of rows
 * @param n Number of columns
 */
void free_mpfr_matrix(mpfr_t ** matrix, int m, int n)
{
    for(int i=0; i < m; i++)
    {
        for(int j=0; j < n; j++)
        {
            mpfr_clear(matrix[i][j]);
        }
        free(matrix[i]);
    }
    free(matrix);
}


/**
 * @brief Transpose a matrix of mpfr_t. The resulting matrix is allocated before calling the function.
 * 
 * @param M Matrix to transpose
 * @param m Number of rows
 * @param n Number of columns
 * @param Mt Resulting matrix
 */
void transpose_mpfr_matrix(mpfr_t **M, int m, int n, mpfr_t ***Mt)
{
    for(int i=0; i < m; i++)
    {
        for(int j=0; j < n; j++)
        {
            mpfr_set((*Mt)[j][i], M[i][j], MPFR_RNDN);
        }
    }
}

/**
 * @brief Copy a matrix of mpfr_t. The resulting matrix is allocated before calling the function.
 * 
 * @param A Matrix to copy
 * @param m Dimension of the matrix
 * @param n Dimension of the matrix
 * @param B Resulting matrix 
 */
void copy_mpfr_matrix(mpfr_t **A, int m, int n, mpfr_t ***B)
{
    for(int i=0; i < m; i++)
    {
        for(int j=0; j < n; j++)
        {
            mpfr_set((*B)[i][j], A[i][j], MPFR_RNDN);
        }
    }
}


void Givens_mpfr(mpfr_t ** A, int m, int n, mpfr_t*** Q, mpfr_t*** R)
{
    // Initialisation of G, Gt, R_tmp, Q_tmp
    mpfr_t **G, **Gt, **R_tmp, **Q_tmp;
    G = (mpfr_t **)malloc(m * sizeof(mpfr_t *));
    Gt = (mpfr_t **)malloc(m * sizeof(mpfr_t *));
    R_tmp = (mpfr_t **)malloc(m * sizeof(mpfr_t *));
    Q_tmp = (mpfr_t **)malloc(m * sizeof(mpfr_t *));
    for(int i=0; i < m; i++)
    {
        G[i] = (mpfr_t *)malloc(m * sizeof(mpfr_t));
        Gt[i] = (mpfr_t *)malloc(m * sizeof(mpfr_t));
        R_tmp[i] = (mpfr_t *)malloc(n * sizeof(mpfr_t));
        Q_tmp[i] = (mpfr_t *)malloc(m * sizeof(mpfr_t));
        for(int j=0; j < m; j++)
        {
            mpfr_init2(G[i][j], 128);
            mpfr_init2(Gt[i][j], 128);
            mpfr_init2(Q_tmp[i][j], 128);
        }
        for(int j=0; j < n; j++)
        {
            mpfr_init2(R_tmp[i][j], 128);
        }
    }

    // Copying A into R
    copy_mpfr_matrix(A, m, n, R);

    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < m; ++i)
        {
            generate_givens_matrix(&G, m, i, j, *R);
            multiply_mpfr_matrices(G, *R, m, m, &R_tmp);
            transpose_mpfr_matrix(G, m, m, &Gt);
            multiply_mpfr_matrices(*Q, Gt, m, m, &Q_tmp);

            // Emptying G
            for(int i=0; i < m; i++)
            {
                for(int j=0; j < m; j++)
                {
                    mpfr_set_d(G[i][j], 0, MPFR_RNDN);
                }
            }

            // Copying R_tmp into R without using copy_matrix_mpfr
            for(int i=0; i < m; i++)
            {
                for(int j=0; j < n; j++)
                {
                    mpfr_set((*R)[i][j], R_tmp[i][j], MPFR_RNDN);
                }
            }

            // Copying Q_tmp into Q without using copy_matrix_mpfr
            for(int i=0; i < m; i++)
            {
                for(int j=0; j < m; j++)
                {
                    mpfr_set((*Q)[i][j], Q_tmp[i][j], MPFR_RNDN);
                }
            }
        }
        
    }

    // Freeing memory
    free_mpfr_matrix(G, m, m);
    free_mpfr_matrix(Gt, m, m);
    free_mpfr_matrix(R_tmp, m, n);
    free_mpfr_matrix(Q_tmp, m, m);
    
}


/**
 * @brief Check if the subdiagonal of a matrix is 0
 * 
 * @param A Matrix to check
 * @param m Number of rows
 * @param n Number of columns
 * @param threshold Precision threshold
 * @return int : 0 if the subdiagonal is 0, 1 otherwise
 */
int thresh_mpfr(mpfr_t ** A, int m, int n, mpfr_t threshold)
{
    mpfr_t tmp;
    mpfr_init2(tmp, 128);
    for(int i=0; i < m; i++)
    {
        for(int j=0; j < n; j++)
        {
            mpfr_abs(tmp, A[i][j], MPFR_RNDN);
            if (mpfr_cmp(tmp, threshold) > 0)
            {
                return 1;
            }
        }
    }
    return 0;
}


/**
 * @brief Recursively compute the eigenvalues of a matrix using the QR algorithm
 * 
 * @param A Matrix to compute the eigenvalues of
 * @param m Number of rows
 * @param n Number of columns
 * @param nb_iter Number of iterations it took to compute the eigenvalues
 * @param threshold Precision threshold
 */
void compute_eigenvalues_mpfr(mpfr_t*** A, int m, int n, int* nb_iter, mpfr_t threshold)
{
    if((thresh_mpfr(*A, m, n, threshold)) == 0)
    {
        return;
    }
    
    // Initializing Q, R, Qt, A2, A_tmp
    mpfr_t **Q, **R, **Qt, **A2, **A_tmp;
    Q       = (mpfr_t **)malloc(m * sizeof(mpfr_t *));
    R       = (mpfr_t **)malloc(m * sizeof(mpfr_t *));
    Qt      = (mpfr_t **)malloc(m * sizeof(mpfr_t *));
    A2      = (mpfr_t **)malloc(m * sizeof(mpfr_t *));
    A_tmp   = (mpfr_t **)malloc(m * sizeof(mpfr_t *));
    for(int i=0; i < m; i++)
    {
        Q[i]        = (mpfr_t *)malloc(m * sizeof(mpfr_t));
        R[i]        = (mpfr_t *)malloc(n * sizeof(mpfr_t));
        Qt[i]       = (mpfr_t *)malloc(m * sizeof(mpfr_t));
        A2[i]       = (mpfr_t *)malloc(n * sizeof(mpfr_t));
        A_tmp[i]    = (mpfr_t *)malloc(n * sizeof(mpfr_t));
        for(int j=0; j < m; j++)
        {
            mpfr_init2(Q[i][j], 128);
            mpfr_init2(Qt[i][j], 128);
        }
        for(int j=0; j < n; j++)
        {
            mpfr_init2(R[i][j], 128);
            mpfr_init2(A2[i][j], 128);
            mpfr_init2(A_tmp[i][j], 128);
        }
    }

    Givens_mpfr(*A, m, n, &Q, &R);
    transpose_mpfr_matrix(Q, m, m, &Qt);
    multiply_mpfr_matrices(Q, *A, m, m, &A2);
    multiply_mpfr_matrices(A2, Qt, m, m, &A_tmp);

    free_mpfr_matrix(Q, m, m);
    free_mpfr_matrix(Qt, m, m);
    free_mpfr_matrix(R, m, n);
    free_mpfr_matrix(A2, m, n);

    (*nb_iter)++;

    for(int i=0; i < m; i++)
    {
        for(int j=0; j < n; j++)
        {
            mpfr_set((*A)[i][j], A_tmp[i][j], MPFR_RNDN);
        }
    }

    free_mpfr_matrix(A_tmp, m, n);

    compute_eigenvalues_mpfr(A, m, n, nb_iter, threshold);
}