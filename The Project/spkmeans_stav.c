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
/* Given N data points (and thier dimension), update the given adj_mat to be the corresponding weighted adjacency matrix.
  if an error occured reteurns FAIL, else- SUCCSESS*/
double** adjacency_matrix(double **data_points, int dimension, int N)
{
    int i, j;
    double **adj_mat = matrix_allocation(N, N);
    if (adj_mat == NULL)
        return NULL;
    for (i = 0; i < N; i++)
    {
        for (j = i; j < N; j++)
        {
            adj_mat[i][j] = (i == j) ? 0 : exp((calc_euclidean_norm(data_points[i], data_points[j], dimension)) / (-2));
            adj_mat[j][i]=adj_mat[i][j];
        }
    }
    return adj_mat;
}

/* Given a vector x, vector y- calculate: ||x-y||2*/
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
    int i, j;
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

    return mat;
}

/* from Zohar, changed names to be more general ask him if OKAY, also fixed error because no c[i][j]=0 was made*/
double **calc_mul(int N, double **A, double **B)
{ 
    /* C=A*B*/
    int i, j, k;
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
    int i, j;
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
    int count;
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
    int i,j;
    i=0;
    j=0;
    /*for(i=0;i<num_rows;i++)
    {
        for (j = 0; j < num_cols; j++)
        {
            if (fscanf(ifp, "%lf", &curr_value) == 1)
            {
                data_input[i][j] = (double)curr_value;
                fgetc(ifp);
            }
            printf("why? %f, \n",curr_value);
        }
    }*/

    /* TODO: this one works for csv too but the other doesnt!, check if it's good..*/
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

/*FOR US...*/
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
                if(mat[i][j]<0 && fabs(mat[i][j]*1000)<1)
                    printf("0.0000");
                else
                    printf("%.4f", mat[i][j]);
            }
            else
            {
                if(mat[i][j]<0 && fabs(mat[i][j])*1000<1)
                    printf("0.0000,");
                else
                    printf("%.4f,", mat[i][j]);
            }
        }
        printf("\n");
    }

}


/* Gets an array to be free (pointers of pointers) and their size*/
void free_memory(double **ArrayToFree, int num_rows){
    int i;
    for (i = 0; i < num_rows; i++){
        free(ArrayToFree[i]);
    }
    free(ArrayToFree);
}


void msg_and_exit(int error_type,int is_error)
{
    if(is_error)
    {
        /* todo handle error- incase of malloc faild or somthing need to exit somehow with error*/
        if (error_type==INVALID_TYPE)
        {
            printf(INVALID);
            exit(1);
        }
        else
        {
            printf(ERROR);
            exit(1);
        }
    }
}

void print_result(double ** mat,int num_rows,int num_cols,enum Goal goal)
{/* TODO add special case for jaacobi?*/
    int i,j;
    /* TODO- added just to avoid errors in moba!!*/
    if(goal==JACOBI)
        num_rows++;
    for (i = 0; i < num_rows; i++)
    {
        for (j = 0; j < num_cols; j++)
        {
            if (j == num_rows - 1)
            {
                if(mat[i][j]<0 && fabs(mat[i][j]*1000)<1)
                    printf("0.0000");
                else
                    printf("%.4f", mat[i][j]);
            }
            else
            {
                if(mat[i][j]<0 && fabs(mat[i][j])*1000<1)
                    printf("0.0000,");
                else
                    printf("%.4f,", mat[i][j]);
            }
        }
        printf("\n");
    }
}


double **run_goal(enum Goal goal,double **data_input,int N,int D)
{
    double **data_output,**wam_matrix,**ddg_matrix;

    /* run WAM*/
    data_output=adjacency_matrix(data_input,D,N);
    if(data_output==NULL)
    {/* an error occured- no need to free*/
        return NULL;
    }
    if(goal==WAM)
        return data_output;
    
    wam_matrix=data_output;
    /* run DDG*/
    data_output=diagonal_matrix(wam_matrix,N);
    if(data_output==NULL)
    {/* an error occured- need to free*/
        free_memory(wam_matrix,N);
        return NULL;
    }
    if(goal==DDG)
        return data_output;

    ddg_matrix=data_output;
    /* run LNORM*/
    data_output=laplacian_matrix(ddg_matrix,wam_matrix,N);
    if(data_output==NULL)
    {/* an error occured- need to free*/
        free_memory(wam_matrix,N);
        free_memory(ddg_matrix,N);
        return NULL;
    }
    if(goal==LNORM)
        return data_output;

    /* TODO for SPK*/
    return NULL;
}

int main(int argc, char *argv[])
{
    char *file_name;
    int N,D;
    double **data_input, **data_output;
    FILE *ifp;
    enum Goal goal=0;

    /* invalid number of arguments*/
    msg_and_exit(INVALID_TYPE,argc!=3);

    /* set goal correct enum*/
    if(!strcmp("wam",argv[1]))
        goal=wam_g;
    if(!strcmp("ddg",argv[1]))
        goal=ddg_g;
    if(!strcmp("lnorm",argv[1]))
        goal=lnorm_g;
    if(!strcmp("jacobi",argv[1]))
        goal=jacobi_g;
    if(!strcmp("spk",argv[1]))
        goal=spk_g;
    msg_and_exit(INVALID_TYPE,goal==0);
    
    file_name=argv[2];
    ifp=fopen(file_name,"r");
    msg_and_exit(ERROR_TYPE,ifp ==NULL);
    
    N=find_N_D(ifp,FIND_N);
    D=find_N_D(ifp,FIND_D);

    /* create matrix for input*/
    data_input=matrix_allocation(N,D);
    msg_and_exit(ERROR_TYPE,data_input==NULL);

    /* set the N points/symetrix matrix in data_input*/
    set_input(ifp,data_input,N,D);
    
    /* set the goal's result in data_output*/
    data_output=run_goal(goal,data_input,N,D);
    if(data_output==NULL)
    { /* an error occured*/
        free_memory(data_input,N);
        msg_and_exit(ERROR_TYPE,1);
    }

    print_result(data_output,N,N,goal);

    fclose(ifp);
    exit(0);
}
