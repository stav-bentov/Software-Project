#ifndef SOFTWARE_PROJECT_SPKMEANS_H
#define SOFTWARE_PROJECT_SPKMEANS_H

#define PY_SSIZE_T_CLEAN
/*#include <Python/Python.h>*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <float.h>

#define ERROR "An Error Has Occurred"
#define INVALID "Invalid Input!"
#define SUCCESS 0
#define FAIL -1
#define FIND_N 1
#define FIND_D 2
#define BUFFER_SIZE 1024
#define INVALID_TYPE 0
#define ERROR_TYPE 1
#define EPSILON 1.0 * pow(10, -5)
#define MAX_ITER 100
#define WAM 1
#define DDG 2
#define LNORM 3
#define JACOBI 4
#define SPK 5
#define INVALID_TYPE 0
#define ERROR_TYPE 1

enum Goal
{ /* TODO: check if SPK is needed*/
  wam_g = 1,
  ddg_g = 2,
  lnorm_g = 3,
  jacobi_g = 4,
  spk_g = 5
};

/* Function's declaretions*/

/* for wam,ddg,lnorm*/
double **adjacency_matrix(double **data_points, int dimension, int N);
double calc_euclidean_norm(double *x, double *y, int dimension);
double **diagonal_matrix(double **adj_mat, int N);
double **laplacian_matrix(double **diag_mat, double **adj_mat, int N);
void cal_D12(double **diag_mat, int N);
double **matrix_allocation(int num_rows, int num_cols);
void check_error(int bool_error);
double **calc_mul(int N, double **A, double **B);
void calc_sub(int N, double **A, double **B);
double **I_matrix(int N);
double **spk_algo(double **lnorm, int N, int K);
void swap(double **mat, int index_1, int index_2);
void sort_matrix_values(double **mat, int l, int r);
int sort_by_p(double **mat, int l, int r);
int eigengap_heuristic(double *eigenvalues, int N);
double **set_T(double **U, int N, int K);

/* for main function*/
double **run_goal(enum Goal goal, double **data_input, int N, int D, int K);
void print_result(double **mat, int num_rows, int num_cols, enum Goal goal);
void msg_and_exit(int error_type, int is_error);
int find_N_D(FILE *ifp, int find_who);
void set_input(FILE *ifp, double **data_input, int num_rows, int num_cols);
void free_memory(double **ArrayToFree, int num_rows);

double **jacobi_algo(int N, double **A);
void matrix_copy(int N, double **dest_mat, double **src_mat);
int check_convergence(int N, double **A, double **A1);
void find_Aij(int N, double **A, int *iPointer, int *jPointer);
void find_c_s_t(double **A, int i, int j, double *cPointer, double *sPointer);
void calc_curr_P(int N, double **curr_P, int i, int j, double c, double s);
void matrix_multiplication(int N, double **src1, double **src2, double **dst);
void get_eigenvalues_from_A1(double *eigenvalues, int N, double **A1);
void transpose(double **mat, int N);
double **jacobi_eigen_merge(int N, double *eigenValues, double **eigenVectors);

#endif