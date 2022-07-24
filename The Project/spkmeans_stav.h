#ifndef SOFTWARE_PROJECT_SPKMEANS_H
#define SOFTWARE_PROJECT_SPKMEANS_H


#define ERROR "An Error Has Occurred"
#define INVALID "Invalid Input!"
#define SUCCESS 0
#define FAIL -1
#define FIND_N 1
#define FIND_D 2
#define BUFFER_SIZE 1024

/* Function's declaretions*/
double **adjacency_matrix(double **data_points,int dimension, int N);
double calc_euclidean_norm(double *x, double *y, int dimension);
double **diagonal_matrix(double **adj_mat, int N);
double **laplacian_matrix(double **diag_mat, double **adj_mat,int N);
void cal_D12(double **diag_mat, int N);
double **matrix_allocation(int num_rows, int num_cols);
void check_error(int bool_error);
double **calc_mul(int N, double **A, double **B);
void calc_sub(int N, double **A,double **B);
double **I_matrix(int N);



int find_N_D( FILE * ifp,int find_who);

#endif