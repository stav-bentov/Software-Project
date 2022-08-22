#ifndef SOFTWARE_PROJECT_SPKMEANS_H
#define SOFTWARE_PROJECT_SPKMEANS_H

#define PY_SSIZE_T_CLEAN
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
#define INVALID_TYPE 0
#define ERROR_TYPE 1
#define EPSILON_JACOBI 0.00001
#define MAX_ITER_JACOBI 100
#define EPSILON_KMEANS 0
#define MAX_ITER_KMEANS 300
#define WAM 1
#define DDG 2
#define LNORM 3
#define JACOBI 4
#define SPK 5
#define SPK_EX2 6

enum Goal
{
  wam_g = 1,
  ddg_g = 2,
  lnorm_g = 3,
  jacobi_g = 4,
  spk_g = 5,
  spk_g2 = 6
};

/* Function's declarations*/
/*Wam,Ddg and Lnorm's functions*/
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
double **spk_algo(double **lnorm, int N, int *K);
double **sort_matrix_values(double **mat, int N);
void eigengap_heuristic(double *eigenvalues, int N, int *K);
double **set_T(double **U, int N, int K);

/* (C) main's functions*/
double **run_goal(enum Goal goal, double **data_input, int N, int D, int *K);
void print_result(double **mat, int num_rows, int num_cols, enum Goal goal);
void msg_and_exit(int error_type, int is_error);
int find_N_D(FILE *ifp, int find_who);
void set_input(FILE *ifp, double **data_input, int num_rows, int num_cols);
void free_memory(double **ArrayToFree, int num_rows);

/* Jacobi's functions*/
double **jacobi_algo(int N, double **A);
void matrix_copy(int num_rows, int num_cols, double **dest_mat, double **src_mat);
int check_convergence(int N, double **A, double **A1);
void find_Aij(int N, double **A, int *iPointer, int *jPointer);
void find_c_s_t(double **A, int i, int j, double *cPointer, double *sPointer);
void calc_curr_P(int N, double **curr_P, int i, int j, double c, double s);
double **jacobi_A_multiplication(int N, double **A, double **A1, double **curr_P, int turn);
void get_eigenvalues_from_A1(double *eigenvalues, int N, double **A1);
void transpose(double **mat, int N);
double **jacobi_eigen_merge(int N, double *eigenValues, double **eigenVectors);

/*Stav jacobi */
void calc_A1(int N, double **A, double c, double s, int i, int j, int *return_value);
double calc_off_diag(int N, double **A);
double **jacobi_algo2(int N, double **A);

/* Kmeans algorithm's functions from ex2*/
int kMeans(int N, int K, double **Datapoints, double **Centroids, int dimension);
int check_euclidean_norm(double **newCentroids, double **oldCentroids, int dimension, int K);
int find_cluster(double **Centroids, double *Datapoint, int dimension, int K);
void updateOldCentroid(double **newCentroids, double **oldCentroids, int dimension, int K);
#endif
