#define PY_SSIZE_T_CLEAN
/*#include <Python/Python.h>*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include "spkmeans_stav.h"

/* 1.1.1 The Weighted Adjacency Matrix*/
/* given N data points (and thier dimension), update the given adj_mat to be the corresponding weighted adjacency matrix.
  if an error occured reteurns FAIL, else- SUCCSESS*/
double** adjacency_matrix(double **data_points, int dimension, int N)
{
    printf("in adjacency_matrix");
    int i, j;
    double **adj_mat = matrix_allocation(N, N);
    if (adj_mat == NULL)
        return NULL;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            adj_mat[i][j] = (i == j) ? 0 : exp((calc_euclidean_norm(data_points[i], data_points[j], dimension)) / (-2));
        }
    }
    return adj_mat;
}

/* given a vector x, vector y- calculate: ||x-y||2*/
double calc_euclidean_norm(double *x, double *y, int dimension)
{
    int j;
    double sum=0;
    for (j = 0; j < dimension; j++)
    {
        sum += pow(x[j] - y[j], 2);
    }
    sum = sqrt(sum);
    return sum;
}

/* 1.1.2 The Diagonal Degree Matrix*/
/* given weighted adjacency matrix size N*N, update the given diag_mat to be the corresponding diagonal degree matrix.
if an error occured reteurns FAIL, else- SUCCSESS*/
double **diagonal_matrix(double **adj_mat, int N)
{
    int i, j, return_value = 0;
    double sum=0;
    double **diag_mat = matrix_allocation(N, N);
    if (diag_mat == NULL)
        return NULL;
    for (i = 0; i < N; i++)
    {
        sum = 0;
        for (j = 0; j < N; j++)
        {
            sum += adj_mat[i][j];
            diag_mat[i][j] = 0;
        }
        diag_mat[i][i] = sum;
    }
    return diag_mat;
}

/* 1.1.3 The Normalized Graph Laplacian*/
double** laplacian_matrix(double **diag_mat, double **adj_mat, int N)
{
    int return_value;
    double **mul1, **mul2;
    double **lnorm;

    cal_D12(diag_mat, N);
    mul1 = calc_mul(N, diag_mat, adj_mat);
    if (mul1 == NULL)
        return NULL;
    mul2 = calc_mul(N, mul1, diag_mat);
    if (mul2 == NULL)
        return NULL;
    lnorm = I_matrix(N);
    if (lnorm == NULL)
        return NULL;

    calc_sub(N, lnorm, mul2);
    return lnorm;
}

/* calculate D^-1/2*/
void cal_D12(double **diag_mat, int N)
{
    int i;
    for (i = 0; i < N; i++)
    {
        diag_mat[i][i] = 1 / sqrt(diag_mat[i][i]);
    }
}

double** matrix_allocation(int num_rows, int num_cols)
{
    /*allocation of memory of size (nxn) and return 1 if there was a failure!*/
    int i;
    double **mat = malloc((sizeof(double *)) * num_rows);
    if (mat == NULL)
        return NULL;
    for (i = 0; i < num_rows; i++)
    {
        mat[i] = malloc((sizeof(double)) * (num_cols));
        if (mat[i] == NULL)
            return NULL;
    }
    
    printf("pointer in matrix_allocation: %p\n",mat);

    return mat;
}

void check_error(int bool_error)
{
    /* todo handle error- incase of malloc faild or somthing need to exit somehow with error*/
    if (bool_error)
    {
        printf("handle error");
    }
}

/* from Zohar, changed names to be more general ask him if OKAY, also fixed error because no c[i][j]=0 was made*/
double **calc_mul(int N, double **A, double **B)
{ 
    /* C=A*B*/
    int i, j, k, return_value;
    double **C = matrix_allocation(N, N);
    if (C == NULL)
        return NULL;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            C[i][j] = 0;
            for (k = 0; k < N; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

/* make subtraction from A and return *updated* A*/
void calc_sub(int N, double **A, double **B)
{ /* A-=B*/
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            A[i][j] -= B[i][j];
        }
    }
}

double **I_matrix(int N)
{
    int i, j, return_value;
    double **I = matrix_allocation(N, N);
    if (I == NULL)
        return NULL;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            I[i][j] = (i == j) ? 1 : 0;
        }
    }
    return I;
}

int find_N_D( FILE * ifp,int find_who)
{/* TODO- under the asumption that in jaacobi also the values sepereted by comas..*/
    int count,first_line;
    char c;

    count=0;

    while ((c = fgetc(ifp)) != EOF)
    {
        if(find_who==FIND_N)
        {
            if (c == '\n')
            {
                count++;
            }
        }
        else
        {
            if (c == '\n')
            {
                rewind(ifp);
                count++;
                return count;
            }
            else
            {
                if (c == ',')
                    count++;
            }
        }
    }

    rewind(ifp);

    if(find_who==FIND_D)
        count++;
    return count;
}

void set_input(FILE *ifp,double **data_input,int num_rows,int num_cols)
{
    printf("in set_input. num_rows: %d, num_cols: %d\n",num_rows,num_cols);
    int i,j;
    i=0;
    j=0;
    /*for(i=0;i<num_rows;i++)
    {
        for (j = 0; j < num_cols; j++)
        {
            //printf("in loop 1");
            if (fscanf(ifp, "%lf", &curr_value) == 1)
            {
                data_input[i][j] = (double)curr_value;
                fgetc(ifp);
            }
            printf("why? %f, \n",curr_value);
        }
    }*/

    // TODO: this one works for csv too but the other doesnt!, check if it's good..
    char buffer[BUFFER_SIZE];
    while (fgets(buffer,BUFFER_SIZE, ifp)) 
    {
        char* curr_value = strtok(buffer, ",");
        while (curr_value) 
        {
            if(j==num_cols)
            {
                j=0;
                i++;
                if(i==num_rows)
                    break;
            }
            data_input[i][j] = atof(curr_value);
            j++;
            curr_value = strtok(NULL, ",");
        }
    }
}

void print_matrix(double ** mat,int num_rows,int num_cols)
{
    int i,j;
    /* writes k-means output file*/
    for (i = 0; i < num_rows; i++)
    {
        for (j = 0; j < num_cols; j++)
        {
            if (j == num_rows - 1)
            {
                if(mat[i][j]<0 && mat[i][j]>-0.00009)
                    printf("0.0000");
                else
                    printf("%.4f", mat[i][j]);
            }
            else
            {
                if(mat[i][j]<0 && mat[i][j]>-0.00009)
                    printf("0.0000,");
                else
                    printf("%.4f,", mat[i][j]);
            }
        }
        printf("\n");
    }

}

#define WAM 1
#define DDG 2
#define LNORM 3
#define JACOBI 4
#define SPK 5

enum goal
{
    wam,
    ddg,
    lnorm,
    jacobi
};

/*static const char *GOAL_ARRAY[]={"wam,ddg,lnorm,jacobi,spk"};*/

/* FOR ME */

int main(int argc, char *argv[])
{
    char *user_goal,*file_name;
    int return_value,N,dimension;
    double **data_input;
    FILE *ifp;

    /* invalid number of arguments*/
    if (argc!=3 || (!strcmp("wam",argv[1]) && !strcmp("ddg",argv[1]) && !strcmp("lnorm",argv[1]) && !strcmp("jacobi",argv[1]) && !strcmp("spk",argv[1])))
    {
        printf(INVALID);
        exit(1);
    }

    user_goal=argv[1];
    file_name=argv[2];

    ifp=fopen(file_name,"r");
    if (ifp ==NULL)
    {
        printf(INVALID);
        exit(1);
    }

    N=find_N_D(ifp,FIND_N);
    dimension=find_N_D(ifp,FIND_D);

    printf("N: %d, D: %d\n",N,dimension);
    printf("%d, %s\n",return_value,user_goal);

    data_input=matrix_allocation(N,dimension);
    if(data_input==NULL)
    {
        printf(ERROR);
        exit(1);
    }
    set_input(ifp,data_input,N,dimension);
    print_matrix(data_input,N,dimension);
    
    double **adj_mat;
    adj_mat=adjacency_matrix(data_input,dimension,N);
    if(adj_mat==NULL)
    {
        printf(ERROR);
        exit(1);
    }
    print_matrix(adj_mat,N,N);

    
    double **diag_mat;
    diag_mat=diagonal_matrix(adj_mat,N);
    if(diag_mat==NULL)
    {
        printf(ERROR);
        exit(1);
    }
    print_matrix(diag_mat,N,N);

    double **lnorm_mat;
    lnorm_mat=laplacian_matrix(diag_mat,adj_mat,N);
    if(lnorm_mat==NULL)
    {
        printf(ERROR);
        exit(1);
    }
    print_matrix(lnorm_mat,N,N);

 
    fclose(ifp);
    return 0;
    exit(0);
}
