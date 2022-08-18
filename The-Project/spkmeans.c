#include "spkmeans.h"

/* ================================== SPK algorithm steps 1-5 ==================================*/
/* Receives N- number of rows, *K- pointer to K (number of clusters) and lnorm's matrix (created from wam,ddg->lnorm)
 * Returns after 3-5 steps from "The Normalized Spectral Clustering Algorithm" a matrix T as N datapoints and calculate K if needed (case K=0)*/
double **spk_algo(double **lnorm, int N, int *K)
{ /* Called after steps 1-2 have been made*/
    double **jacobi_output, **U, **eigenvectors, **T, **sort_transpose;

    jacobi_output = (double **)jacobi_algo(N, lnorm);
    if (jacobi_output == NULL)
        return NULL;

    /* Transpose on eigenvectors- to make the sort easier*/
    eigenvectors = jacobi_output + 1; /* jacobi without eigenvalues*/
    transpose(eigenvectors, N);
    /*sort_matrix_values(jacobi_output, N);*/
    sort_transpose= sort_matrix_values(jacobi_output,N);
    if(sort_transpose == NULL)
        return NULL;

    /* The Eigengap Heuristic- was told not to handle a case where k=1*/
    if (*K == 0)
        eigengap_heuristic(sort_transpose[0], N, K);

    eigenvectors = sort_transpose + 1; /* jacobi without eigenvalues*/
    transpose(eigenvectors, N);

    /* U points to the start of eigenvectors, we will use only the first K vectors (first K columns)*/
    U = eigenvectors;
    T = set_T(U, N, *K);
    if (T == NULL)
    {
        free_memory(jacobi_output, N + 1);
        return NULL;
    }
    return T;
}

double** sort_matrix_values(double **mat, int N)
{
    int i, j, max_index;
    double max_value;
    double **sort_mat =matrix_allocation(N+1,N);
    if(sort_mat == NULL)
    {
        free_memory(mat,N+1);
        return NULL;
    }
    for (i = 0; i < N; i++)
    {
        max_index = -1;
        max_value =-1;
        for (j = 0; j < N; j++)
        {
            /* found new max*/
            if (max_value < mat[0][j])
            {
                max_index = j;
                max_value=mat[0][j];
            }
        }
        sort_mat[0][i]=max_value;
        sort_mat[i+1]=mat[max_index+1];
        mat[0][max_index]=-1;
    }
    return sort_mat;
}

/* Receives a jacobi's matrix
 * Sort first row and rows 1 to N according to the eigenvalues in first row */
/*void sort_matrix_values(double **mat, int N, int *return_value)
{
    int i, j, max_index,swap_index,temp_index;
    double max_value,temp_value;
    double *sort_values =calloc(N,sizeof(double));
    if(sort_values == NULL)
    {
        *return_value=0;
        return;
    }
    double *sort_vectors =malloc(N*sizeof(*double));
    if(sort_vectors == NULL)
    {
        *return_value=0;
        free(sort_values);
        return;
    }
    for (i = 0; i < N; i++)
    {
        max_index = -1;
        max_value = -1;
        for (j = 0; j < N; j++)
        {
            if (max_value < mat[0][j])
            {
                max_index = j;
                max_value=mat[0][j];
            }
        }
        sort_values[i]=max_value;
        sort_vectors[i]=mat[max_index];
        mat[0][max_index]=-1;
    }

    for (i=1; i<N+1; i++)
    {
        mat[0][i-1]=sort_values[i-1];
        mat[i]=sort_vectors[i-1];
    }
    free(sort_values);
    free(sort_vectors);
}*/

/* Receives a jacobi's matrix
 * Sort first row and rows 1 to N according to the eigenvalues in first row */
/*void sort_matrix_values(double **mat, int N)
{
    int i, j, i_max_swap;
    double max_value;
    for (i = 0; i < N; i++)
    {
        i_max_swap = -1;
        max_value = mat[0][i];
        for (j = i + 1; j < N; j++)
        {
            if (max_value < mat[0][j])
            {
                i_max_swap = j;
                max_value = mat[0][j];
            }
        }
        if (i_max_swap != -1)
            swap(mat, i_max_swap, i);
    }
}*/

/* Receives jacobi's matrix, index_1 and index_2 between 0 to N-1
 * Swaps between mat[0][index_1] to mat[0][index_2] and mat[index_1+1] to mat[index_2+1]*/
void swap(double **mat, int index_1, int index_2)
{
    double temp_value, *temp_vector;

    temp_value = mat[0][index_1];
    mat[0][index_1] = mat[0][index_2];
    mat[0][index_2] = temp_value;

    temp_vector = mat[index_1 + 1];
    mat[index_1 + 1] = mat[index_2 + 1];
    mat[index_2 + 1] = temp_vector;
}

/* Receives U (created by largest eigenvectors of jacobi), N- number of rows, K- number of columns
 * Returns T- by renormalizing each of U’s rows to have unit length */
double **set_T(double **U, int N, int K)
{
    int i, j, q;
    double sum=0;

    double **T = matrix_allocation(N, K);
    if (T == NULL)
        return NULL;
    
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K; j++)
        {
            if (j == 0)
            { /* Calculate sum once for each new row!*/
                sum = 0;
                for (q = 0; q < K; q++)
                {
                    sum += pow(U[i][q], 2);
                }
            }
            if(sum != 0)
                T[i][j] = U[i][j] / sqrt(sum);
            else
                return NULL;
        }
    }
    return T;
}

/* Receives eigenvalues (in decreasing order), N- number of values, K- number of clusters
 * Calculate K as described (Eigengap Heuristic algorithm) */
void eigengap_heuristic(double *eigenvalues, int N, int *K)
{ /* (Assumption) lnorm formed as a decreasing ordered eigenvalues*/
    int i;
    double curr_max_gap = DBL_MIN;

    /* lmda(1)= E[0]>=lmda(2)=E[1]>=...>=lmda(n/2)=E[(N/2)-1]>=0*/
    for (i = 1; i <= (int)(N / 2); i++)
    {
        if (curr_max_gap < fabs(eigenvalues[i - 1] - eigenvalues[i]))
        {
            curr_max_gap = fabs(eigenvalues[i - 1] - eigenvalues[i]);
            *K = i;
        }
    }
}
/* ================================== Done SPK ==================================*/

/* ================================== WAM (Weighted Adjacency Matrix) ================================== */
/* Receives N datapoints and their dimension 
 * Returns the corresponding weighted adjacency matrix.
 * If an error occurred returns NULL*/
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
            adj_mat[i][j] = (i == j) ? 0 : (exp((calc_euclidean_norm(data_points[i], data_points[j], dimension)) / (-2)));
            adj_mat[j][i] = adj_mat[i][j];
        }
    }
    return adj_mat;
}

/* Receives 2 vectors- x,y and thier dimension
 * Returns their distance: ||x-y||2*/
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
/* ================================== Done WAM ==================================*/

/* ================================== DDG (Diagonal Degree Matrix) ================================== */
/* Receives weighted adjacency matrix and N- number of rows/columns
 * Returns the corresponding diagonal degree matrix.
 * If an error occurred returns NULL*/
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
/* ================================== Done DDG ==================================*/

/* ================================== LNORM (Normalized Graph Laplacian) ================================== */
/* Receives diagonal degree matrix,weighted adjacency matrix and N- number of rows/columns
 * Returns the corresponding normalized graph Laplacian matrix.
 * If an error occurred returns NULL*/
double **laplacian_matrix(double **diag_mat, double **adj_mat, int N)
{
    /* mul1 = diag_mat^(-0.5)*adj_mat^(-0.5) 
     * mul2 = mul1*diag_mat^(-0.5) = diag_mat^(-0.5)*adj_mat*diag_mat^(-0.5) 
     * lnorm=I - mul2 = I - diag_mat^(-0.5)*adj_mat*diag_mat^(-0.5) */
    double **mul1, **mul2, **lnorm;

    cal_D12(diag_mat, N);
    mul1 = calc_mul(N, diag_mat, adj_mat);
    if (mul1 == NULL)
        return NULL;
    mul2 = calc_mul(N, mul1, diag_mat);
    if (mul2 == NULL)
    {
        free_memory(mul1, N);
        return NULL;
    }
    lnorm = I_matrix(N);
    if (lnorm == NULL)
    {
        free_memory(mul1, N);
        free_memory(mul2, N);
        return NULL;
    }
    calc_sub(N, lnorm, mul2);
    return lnorm;
}

/* Receives D- diagonal degree matrix, N- number of rows/columns
 * Updates D to D^(-0.5)*/ 
void cal_D12(double **diag_mat, int N)
{
    int i;
    for (i = 0; i < N; i++)
    { /* TODO check if Rami answered :https://moodle.tau.ac.il/mod/forum/discuss.php?d=128591*/
        /* (assumption) diag_mat[i][i]!=0 */
        diag_mat[i][i] = 1 / sqrt(diag_mat[i][i]);
    }
}

/* Receives matrixs: A,B and N- number of rows/columns
 * Returns matrix C= A*B
 * If an error occurred returns NULL*/
double **calc_mul(int N, double **A, double **B)
{
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

/* Receives matrixs: A,B and N- number of rows/columns
 * Updates A= A-B */
void calc_sub(int N, double **A, double **B)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            A[i][j] -= B[i][j];
        }
    }
}

/* Receives a number-N
 * Returns the Identity matrix from size N*N */
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
/* ================================== Done LNORM ==================================*/

/* ================================== JACOBI ================================== */
double **jacobi_algo(int N, double **A)
{
    int counter,iPointer, jPointer;
    double cPointer, sPointer;                   /*pivot element, s,c*/
    double **A1, **V, **curr_P, **jacobi_result; /* A' matrix, eigenVectors, *P matrix - keeps changing and (V = V x curr_P)*  */
    double *eigenvalues;

    counter=0;

    /* Memory allocations- if matrix_allocation returns NULL-
     * an error occurred- free previous allocations and return NULL */
    A1 = matrix_allocation(N, N);
    if (A1 == NULL)
        return NULL;

    V = I_matrix(N);
    if (V == NULL)
    {
        free_memory(A1, N);
        return NULL;
    }

    curr_P = matrix_allocation(N, N);
    if (curr_P == NULL)
    {
        free_memory(A1, N);
        free_memory(V, N);
        return NULL;
    }

    eigenvalues = malloc(N * sizeof(double)); /*len of diagonal of squared matrix (NxN) is always N*/
    if (eigenvalues == NULL)
    {
        /*TODO decide which free memory to use 1 or regular*/
        free_memory(A1, N);
        free_memory(V, N);
        free_memory(curr_P, N);
        /*free_memory1(N, 3, A1, V, curr_P);*/
        return NULL;
    }

    while ((MAX_ITER_JACOBI > counter) && (!check_convergence(N, A, A1) || (counter == 0)))
    {
        counter++;
        /*A = A1*/
        if (counter != 1)
            matrix_copy(N, N, A, A1);

        find_Aij(N, A, &iPointer, &jPointer);
        find_c_s_t(A, iPointer, jPointer, &cPointer, &sPointer);
        calc_curr_P(N, curr_P, iPointer, jPointer, cPointer, sPointer);

        A1 = jacobi_A_multiplication(N, A, A1, curr_P, 0);
        if (A1 == NULL)
        {
            free(eigenvalues);
            free_memory(V, N);
        }
        A1 = jacobi_A_multiplication(N, A, A1, curr_P, 1);
        if (A1 == NULL)
        {
            free(eigenvalues);
            free_memory(V, N);
        }

        V = calc_mul(N, V, curr_P);
        if (V == NULL)
        { /* it is reachable!*/
            free(eigenvalues);
            free_memory(A1, N);
            free_memory(curr_P, N);
            /*free_memory1(N, 2, A1, curr_P);*/
            return NULL;
        }
    }

    get_eigenvalues_from_A1(eigenvalues, N, A1); /*Getting eigenvalues from A' diagonal!*/
    jacobi_result = jacobi_eigen_merge(N, eigenvalues, V);
    if (jacobi_result == NULL)
    {
        free_memory(A1, N);
        free_memory(V, N);
        free_memory(curr_P, N);
        /*free_memory1(N, 3, A1, V, curr_P);*/
        free(eigenvalues);
    }

    return jacobi_result;
}

double **jacobi_A_multiplication(int N, double **A, double **A1, double **curr_P, int turn)
{
    double **A1_pointer;
    double **mat1, **mat2;

    mat1 = (turn == 0) ? curr_P : A1;
    mat2 = (turn == 0) ? A : curr_P;

    transpose(curr_P, N);
    A1_pointer = A1;
    A1 = calc_mul(N, mat1, mat2);
    free_memory(A1_pointer, N);
    if (A1 == NULL)
        free_memory(curr_P, N);

    return A1;
}

/* Receives matrix mat and N- number of rows/columns
 * Updates mat to be transpose(mat) */
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
    if (off_A_squared - off_A1_squared <= EPSILON_JACOBI)
        return 1;
    return 0;
}

void find_Aij(int N, double **A, int *iPointer, int *jPointer)
{ /* finds the off-diagonal element with the largest ABSOLUTE value*/
    int q, l;
    double maxElem = -DBL_MAX; /*todo change - be careful*/

    for (q = 0; q < N; ++q)
    {
        for (l = q + 1; l < N; ++l)
        {
            if (fabs(A[q][l]) > maxElem)
            {
                maxElem = fabs(A[q][l]);
                /*changes i,j pointers*/
                *iPointer = q;
                *jPointer = l;
            }
        }
    }
}

void find_c_s_t(double **A, int i, int j, double *cPointer, double *sPointer)
{
    double theta, t;
    double signTheta = 1;
    if (A[i][j] == 0){
        *cPointer = 1;
        *sPointer = 0;
        return;
    }
    theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);

    if (theta < 0)
        signTheta = -1;
    t = (signTheta) / (fabs(theta) + sqrt(theta * theta + 1));
    *cPointer = (1) / (sqrt(1 + t * t));
    *sPointer = t / sqrt(1 + t * t);
}

void calc_curr_P(int N, double **curr_P, int i, int j, double c, double s)
{
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
    int i;
    res = matrix_allocation(N + 1, N);
    if (res == NULL)
    {
        return NULL;
    }

    for (i = 0; i < N; ++i)
        res[0][i] = eigenValues[i];
    matrix_copy(N, N, &res[1], eigenVectors);

    return res;
}
/* ================================== Done JACOBI ==================================*/

/* ================================== General/ Main's Functions ==================================*/
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

/* Receives file's pointer and an int- find_who (2 options: N (FIND_N) or D (FIND_D))
 * Returns N or D according to find_who */
int find_N_D(FILE *ifp, int find_who)
{
    int count;
    char c;

    count = 0;

    while ((c = fgetc(ifp)) != EOF)
    {   
        /* N- calculated by number of rows ("\n")*/
        if (find_who == FIND_N)
        {
            if (c == '\n')
            {
                count++;
            }
        }
        else
        {
            /* D- calculated by number of comas (+ 1) in first row.
             * After done reading first row- returns calculated D */
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

    /*TODO delete this!
    if (find_who == FIND_D)
        count++;*/
        
    return count;
}

/* Receives matrixs dest_mat,src_mat and number of rows and columns
 * Updates dest_mat to be a copy of src_mat */
void matrix_copy(int num_rows, int num_cols, double **dest_mat, double **src_mat)
{
    int i, j;
    for (i = 0; i < num_rows; i++)
    {
        for (j = 0; j < num_cols; j++)
        {
            dest_mat[i][j] = src_mat[i][j];
        }
    }
}

/* Receives file's pointer, pointer to an empty matrix of datapoints and num of rows, num of cols
 * Updated data_input according to the file's data */
void set_input(FILE *ifp, double **data_input, int num_rows, int num_cols)
{
    int i, j;
    double curr_value;
    i = 0;
    j = 0;
    for (i = 0; i < num_rows; i++)
    {
        for (j = 0; j < num_cols; j++)
        {
            if (fscanf(ifp, "%lf", &curr_value) == 1)
                data_input[i][j] = curr_value;
            else
            {
                /*TODO check with weird input files!!!!!!!*/
                j--;
            }
            fgetc(ifp);
        }
    }
}

/* Receives a matrix: mat_to_free and num_rows- matrix's number of rows 
 * Frees mat_to_free*/
void free_memory(double **mat_to_free, int num_rows)
{
    int i;
    for (i = 0; i < num_rows; i++)
    {
        free(mat_to_free[i]);
    }
    free(mat_to_free);
}

/* Receives an int: error_type (2 options: invalid (INVALID_TYPE) or error (ERROR_TYPE)) and is_error- boolean number
 * If is_error is 1 (true- an error occurred), print the correspond message (according to error_type) and exit
 * Else- continue (do nothing) */
void msg_and_exit(int error_type, int is_error)
{
    if (is_error)
    {
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

/* Receives a matrix, number of rows and columns, and enum Goal
 * Prints the matrix (if the goal is jacobi then updates num_rows +1) */
void print_result(double **mat, int num_rows, int num_cols, enum Goal goal)
{
    int i, j;

    if (goal == JACOBI)
        num_rows++;

    for (i = 0; i < num_rows; i++)
    {
        for (j = 0; j < num_cols; j++)
        {
            if (j == num_cols - 1)
            {
                printf("%.4f", mat[i][j]);
            }
            else
            {
                printf("%.4f,", mat[i][j]);
            }
        }
        printf("\n");
    }
}


/* Receives an enum Goal, matrix with given data (from file), N- number of rows (or number of points), 
 * D- number of columns (or point's dimension) and a pointer to K
 * Returns corresponed matrix according to the given goal
 * If an error occured returns NULL*/
double **run_goal(enum Goal goal, double **data_input, int N, int D, int *K)
{
    double **data_output, **wam_matrix, **ddg_matrix, **lnorm_matrix;

    if (goal == JACOBI)
    {
        data_output = jacobi_algo(N, data_input);
        return data_output;
    }

    /* Run WAM*/
    data_output = adjacency_matrix(data_input, D, N);
    if (data_output == NULL)
    { /* An error occurred- no need to free*/
        return NULL;
    }
    if (goal == WAM)
        return data_output;

    wam_matrix = data_output;

    /* Run DDG*/
    data_output = diagonal_matrix(wam_matrix, N);
    if (data_output == NULL)
    { /* An error occurred- need to free*/
        free_memory(wam_matrix, N);
        return NULL;
    }
    if (goal == DDG)
        return data_output;

    ddg_matrix = data_output;

    /* Run LNORM*/
    data_output = laplacian_matrix(ddg_matrix, wam_matrix, N);
    if (data_output == NULL)
    { /* An error occurred- need to free*/
        free_memory(wam_matrix, N);
        free_memory(ddg_matrix, N);
        return NULL;
    }
    if (goal == LNORM)
        return data_output;

    lnorm_matrix = data_output;

    /* run SPK*/
    data_output = spk_algo(lnorm_matrix, N, K);
    if (data_output == NULL)
    { /* An error occurred- need to free*/
        free_memory(wam_matrix, N);
        free_memory(ddg_matrix, N);
        free_memory(lnorm_matrix, N);
        return NULL;
    }
    return data_output;
}

/* Receives k, goal and file_name from user
 * Calculates needed information from file and call run_goal, prints result at end*/
int main(int argc, char *argv[])
{
    char *file_name;
    int N, D, K; /*todo- delete K*/
    double **data_input, **data_output;
    FILE *ifp;
    enum Goal goal = 0;
    K = 0;

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
    msg_and_exit(INVALID_TYPE, goal == 0);

    file_name = argv[2];
    ifp = fopen(file_name, "r");
    msg_and_exit(ERROR_TYPE, ifp == NULL);

    N = find_N_D(ifp, FIND_N);
    D = find_N_D(ifp, FIND_D);

    /* Creates matrix for input*/
    data_input = matrix_allocation(N, D);
    msg_and_exit(ERROR_TYPE, data_input == NULL);

    /* Sets the N points/symmetric matrix in data_input*/
    set_input(ifp, data_input, N, D);

    /* Sets the goal's result in data_output*/
    data_output = run_goal(goal, data_input, N, D, &K);
    if (data_output == NULL)
    { /* An error has occurred*/
        free_memory(data_input, N);
        msg_and_exit(ERROR_TYPE, 1); /* todo check if \n is not necessary after error message */
    }
    
    print_result(data_output, N, N, goal);

    fclose(ifp);
    exit(0);
}
