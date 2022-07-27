#define PY_SSIZE_T_CLEAN
/*#include <Python/Python.h>*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include "spkmeans.h"

/* 1.1.1 The Weighted Adjacency Matrix*/
/* Given N data points (and thier dimension), update the given adj_mat to be the corresponding weighted adjacency matrix.
  if an error occured reteurns FAIL, else- SUCCSESS*/
double **adjacency_matrix(double **data_points, int dimension, int N)
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
            adj_mat[j][i] = adj_mat[i][j];
        }
    }
    return adj_mat;
}

/* Given a vector x, vector y- calculate: ||x-y||2*/
double calc_euclidean_norm(double *x, double *y, int dimension)
{
    int j;
    double sum = 0;
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
    double sum = 0;
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
double **laplacian_matrix(double **diag_mat, double **adj_mat, int N)
{ // TODO: handle memory- match with zohar
    double **mul1, **mul2, **lnorm;

    mul1=matrix_allocation(N,N);
    if(mul1==NULL)
        return NULL;
    mul2=matrix_allocation(N,N);
    if(mul2==NULL)
    {
        free_memory(mul1,N);
        return NULL;
    }
    lnorm = I_matrix(N);
    if (lnorm == NULL)
        return NULL;

    cal_D12(diag_mat, N);
    matrix_multiplication(N,diag_mat,adj_mat,mul1);
    if(mul1==NULL)
        return NULL;
    matrix_multiplication(N,mul1,diag_mat,mul2);
    if(mul2==NULL)
    {
        free_memory(mul2,N);
        return NULL;
    }
    calc_sub(N, lnorm, mul2);
    return lnorm;

    //mul1 = calc_mul(N, diag_mat, adj_mat);
    //mul2 = calc_mul(N, mul1, diag_mat);
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

double **matrix_allocation(int num_rows, int num_cols)
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

double **jacobi_algo(int N, double **A)
{/* TODO: check if can erase NULL*/
    int counter = 0;
    int iPointer, jPointer;
    double Aij;                        /*pivot element*/
    double **A1 = NULL;                /* A' matrix */
    double **V = NULL, **V_tag = NULL; /*eigenVectors*/
    double **curr_P,**jacobi_result;   /*P matrix - keeps changing and (V = V x curr_P)*/
    double *eigenvalues;
    double cPointer, sPointer;

    A1 = matrix_allocation(N, N);
    if (A1 == NULL)
    { /* an error occured- no need to free*/
        return NULL;
    }

    V = I_matrix(N);
    if (V == NULL)
    {
        free_memory(A1, N);
        return NULL;
    }

    curr_P = matrix_allocation(N, N);
    if (curr_P == NULL)
    { /* an error occured- need to free*/
        free_memory(A1, N);
        free_memory(V, N);
        return NULL;
    }

    eigenvalues = malloc(N * sizeof(double)); /*len of diagonal of squared matrix (NxN) is always N*/
    if (eigenvalues == NULL)
    {
        free_memory(A1, N);
        free_memory(V, N);
        free_memory(curr_P, N);
        return NULL;
    }

    while ((MAX_ITER > counter) && (!check_convergence(N, A, A1) || (counter == 0)))
    { /*todo check : even if check Convergence is true in the beginning, i want to do this while loop*/
        counter++;
        /*A = A1*/
        if (counter != 1)
            matrix_copy(N, A, A1);

        Aij = find_Aij(N, A, &iPointer, &jPointer);
        find_c_s_t(A, Aij, iPointer, jPointer, &cPointer, &sPointer);
        calc_curr_P(N, curr_P, iPointer, jPointer, cPointer, sPointer);
        /*calc_A1(N, A, A1,cPointer,sPointer,iPointer,jPointer);gets A, c, s, i, j*/

        /*TODO: change to one func*/
        transpose(curr_P, N);                     /* P -> P_transpose */
        matrix_multiplication(N, curr_P, A, A1);  /* A' = P_transpose*A */
        transpose(curr_P, N);                     /* P_transpose -> P */
        matrix_multiplication(N, A1, curr_P, A1); /* A' = (P_transpose*A)*P */
        matrix_multiplication(N, V, curr_P, V);   /* V = V * curr_P*/
    }

    get_eigenvalues_from_A1(eigenvalues, N, A1); /*getting eigenvalues from A' diagonal!*/
    jacobi_result=jacobi_eigen_merge(N, eigenvalues, V);
    if(jacobi_result==NULL)
    {
        free_memory(A1, N);
        free_memory(V, N);
        free_memory(curr_P, N);
        //free_memory(eigenvalues, N);
        // TODO: check if change eigenvalues allocation changes this!
        free(eigenvalues);
    }
    return jacobi_result;
}

void transpose(double **mat, int N)
{
    int i, j;
    double tmp;

    for (i = 0; i < N; i++)
    {
        for (j = i + 1; j < N; j++)
        {
            tmp = mat[i][j];
            mat[i][j] = mat[j][i];
            mat[j][i] = tmp;
        }
    }
}

void matrix_copy(int N, double **dest_mat, double **src_mat)
{
    int i, j;

    for (i = 0; i < N; i++)
    { /* rad-nr */
        for (j = 0; j < N; j++)
        { /* kolumn-nr */
            dest_mat[i][j] = src_mat[i][j];
        }
    }
}

void matrix_multiplication(int N, double **src1, double **src2, double **dst)
{
    double **temp;
    int i, j, k;
    temp = matrix_allocation(N, N);
    if (temp == NULL)
    {
        free_memory(dst, N);
        dst = NULL;
        return;
    }

    for (i = 0; i < N; i++)
    { /* temp = src1*src2 */
        for (j = 0; j < N; j++)
        {
            temp[i][j] = 0;
            for (k = 0; k < N; k++)
                temp[i][j] += src1[i][k] * src2[k][j];
        }
    }
    for (i = 0; i < N; i++)
    { /* dst = temp */
        for (j = 0; j < N; j++)
        {
            dst[i][j] = temp[i][j];
        }
    }

    free_memory(temp, N);
}

int check_convergence(int N, double **A, double **A1)
{
    /* (true) iff off(A)^2 - off(A')^2 <= epsilon */
    int i, j;
    double off_A_squared = 0, off_A1_squared = 0;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            off_A_squared += i == j ? 0 : pow(A[i][j], 2);
            off_A1_squared += i == j ? 0 : pow(A1[i][j], 2);
        }
    }
    if (off_A_squared - off_A1_squared <= EPSILON)
        return 1;
    return 0;
}

double find_Aij(int N, double **A, int *iPointer, int *jPointer)
{ /* finds the off-diagonal element with the largest ABSOLUTE value*/
    int q, l;
    double maxElem = -DBL_MAX; /*todo change - be careful*/

    for (q = 0; q < N; ++q)
    {
        for (l = 0; l < N; ++l)
        {
            if ((q != l) && fabs(A[q][l]) > maxElem)
            {
                maxElem = fabs(A[q][l]);
                /*changes i,j pointers*/
                *iPointer = q;
                *jPointer = l;
            }
        }
    }
    return maxElem; /*value of Aij*/
}

void find_c_s_t(double **A, double aij, int i, int j, double *cPointer, double *sPointer)
{
    double theta, t;
    double signTheta = 1; /*todo check if sign(theta) is 0 / 1 or something else*/
    theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);

    if (theta < 0)
        signTheta = -1;
    t = (signTheta) / (fabs(theta) + sqrt(theta * theta + 1));
    *cPointer = (1) / (sqrt(1 + t * t));
    *sPointer = t / sqrt(1 + t * t);
}

void calc_curr_P(int N, double **curr_P, int i, int j, double c, double s)
{ /*todo check, maybe this function call is not necessary - complexity wise and no return val*/
    int k, l;

    for (k = 0; k < N; ++k)
    {
        for (l = 0; l < N; ++l)
        {
            if (k == l)
                curr_P[k][l] = (k == i || l == j) ? c : 1;
            else if (k == j && l == i)
                curr_P[k][l] = -s;
            else
                curr_P[k][l] = (k == i && l == j) ? s : 0;
        }
    }
}

void get_eigenvalues_from_A1(double *eigenvalues, int N, double **A1)
{
    int i;
    for (i = 0; i < N; ++i)
        eigenvalues[i] = A1[i][i];
}

double **jacobi_eigen_merge(int N, double *eigenValues, double **eigenVectors)
{
    double **res = NULL;
    int i, j;
    res = matrix_allocation(N + 1, N);
    if(res==NULL)
    {
        return NULL;
    }

    for (i = 0; i < N; ++i)
        res[0][i] = eigenValues[i];

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            res[i + 1][j] = eigenVectors[i][j];
        }
    }
    return res;
}

/* MAIN's functions*/
int find_N_D(FILE *ifp, int find_who)
{ /* TODO- under the asumption that in jaacobi also the values sepereted by comas..*/
    int count;
    char c;

    count = 0;

    while ((c = fgetc(ifp)) != EOF)
    {
        if (find_who == FIND_N)
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

    if (find_who == FIND_D)
        count++;
    return count;
}

void set_input(FILE *ifp, double **data_input, int num_rows, int num_cols)
{
    int i, j;
    i = 0;
    j = 0;
    // TODO: checked- doesnt work on csv
    /*double curr_value;
    for(i=0;i<num_rows;i++)
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
    while (fgets(buffer, BUFFER_SIZE, ifp))
    {
        char *curr_value = strtok(buffer, ",");
        while (curr_value)
        {
            if (j == num_cols)
            {
                j = 0;
                i++;
                if (i == num_rows)
                    break;
            }
            data_input[i][j] = atof(curr_value);
            j++;
            curr_value = strtok(NULL, ",");
        }
    }
}

/* Gets an array to be free (pointers of pointers) and their size*/
void free_memory(double **ArrayToFree, int num_rows)
{
    int i;
    for (i = 0; i < num_rows; i++)
    {
        free(ArrayToFree[i]);
    }
    free(ArrayToFree);
}

void msg_and_exit(int error_type, int is_error)
{
    if (is_error)
    {
        /* todo handle error- incase of malloc faild or somthing need to exit somehow with error*/
        if (error_type == INVALID_TYPE)
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

void print_result(double **mat, int num_rows, int num_cols, enum Goal goal)
{
    int i, j;
    if (goal == JACOBI)
        num_rows++;
    for (i = 0; i < num_rows; i++)
    {
        for (j = 0; j < num_cols; j++)
        {
            if (j == num_rows - 1)
            {
                if (mat[i][j] < 0 && fabs(mat[i][j] * 1000) < 1)
                    printf("0.0000");
                else
                    printf("%.4f", mat[i][j]);
            }
            else
            {
                if (mat[i][j] < 0 && fabs(mat[i][j]) * 1000 < 1)
                    printf("0.0000,");
                else
                    printf("%.4f,", mat[i][j]);
            }
        }
        printf("\n");
    }
}

double **run_goal(enum Goal goal, double **data_input, int N, int D)
{
    double **data_output, **wam_matrix, **ddg_matrix;

    if (goal == JACOBI)
    {
        data_output = jacobi_algo(N,data_input);
        return data_output;
    }

    /* run WAM*/
    data_output = adjacency_matrix(data_input, D, N);
    if (data_output == NULL)
    { /* an error occured- no need to free*/
        return NULL;
    }
    if (goal == WAM)
        return data_output;

    wam_matrix = data_output;
    /* run DDG*/
    data_output = diagonal_matrix(wam_matrix, N);
    if (data_output == NULL)
    { /* an error occured- need to free*/
        free_memory(wam_matrix, N);
        return NULL;
    }
    if (goal == DDG)
        return data_output;

    ddg_matrix = data_output;
    /* run LNORM*/
    data_output = laplacian_matrix(ddg_matrix, wam_matrix, N);
    if (data_output == NULL)
    { /* an error occured- need to free*/
        free_memory(wam_matrix, N);
        free_memory(ddg_matrix, N);
        return NULL;
    }
    if (goal == LNORM)
        return data_output;

    /* TODO for SPK*/
    return NULL;
}

int main(int argc, char *argv[])
{
    char *file_name;
    int N, D;
    double **data_input, **data_output;
    FILE *ifp;
    enum Goal goal = 0;

    /* invalid number of arguments*/
    msg_and_exit(INVALID_TYPE, argc != 3);
    
    /* set goal correct enum*/
    if (!strcmp("wam", argv[1]))
        goal = wam_g;
    if (!strcmp("ddg", argv[1]))
        goal = ddg_g;
    if (!strcmp("lnorm", argv[1]))
        goal = lnorm_g;
    if (!strcmp("jacobi", argv[1]))
        goal = jacobi_g;
    if (!strcmp("spk", argv[1]))
        goal = spk_g;
    msg_and_exit(INVALID_TYPE, goal == 0);

    file_name = argv[2];
    ifp = fopen(file_name, "r");
    msg_and_exit(ERROR_TYPE, ifp == NULL);

    N = find_N_D(ifp, FIND_N);
    D = find_N_D(ifp, FIND_D);

    /* create matrix for input*/
    data_input = matrix_allocation(N, D);
    msg_and_exit(ERROR_TYPE, data_input == NULL);

    /* set the N points/symetrix matrix in data_input*/
    set_input(ifp, data_input, N, D);

    /* set the goal's result in data_output*/
    data_output = run_goal(goal, data_input, N, D);
    if (data_output == NULL)
    { /* an error occured*/
        free_memory(data_input, N);
        msg_and_exit(ERROR_TYPE, 1);
    }

    print_result(data_output, N, N, goal);

    fclose(ifp);
    exit(0);
}
