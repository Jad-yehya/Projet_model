//
// Created by Jad Yehya and Ezzedine Chahine on 24/10/2022.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef MODEL_PROJECT_TOOLS_H
#define MODEL_PROJECT_TOOLS_H

void create_matrix(double **A, int n);
void swap_row(double **A, int i, int j);
void swap_columns(double **A, int i, int j, int n);
void print_matrix(double **A,int m, int n);
double** multiply_matrices(double **A, double **B, int n);
void findGreatest(double **A, int i, int j, int n,  int *in, int *jn);
void generate_G(double **G, int i, int j, int m, double** R);
void transpose_matrix(double ** M, int m, int n, double*** Mt);
void Givens(double **A, int m, int n, double** Q, double **R);
void Givens2(double **A, int m, int n, double*** Q, double*** R);
void copy(double **A, double **B, int m, int n);
void compute_eigenvalues(double ***A, int m, int n, int* nb_iter);
void compute_eigenvalues2(double **A, int m, int n);
void free_matrix(double **A, int m);
void multiplymatrices(double **A, double **B, int n, double ***C);
void compute_eigenvalues3(double ***A, int m, int n, long* nb_iter);
void Givens3(double **A, int m, int n, double*** Q, double*** R);
void generate_random_matrix(double ***A, int m, int n, int a, int b);
void generate_symmetric_matrix(double ***A, int m, int n, int a, int b);

#endif //MODEL_PROJECT_TOOLS_H
