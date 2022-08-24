#include "spkmeans.h"

/* ================================== SPK algorithm steps 1-5 ==================================*/

/* Receives N- number of rows, *K- pointer to K (number of clusters) and lnorm's matrix (created from wam,ddg->lnorm)
 * Returns after 3-5 steps from "The Normalized Spectral Clustering Algorithm" a matrix T as N datapoints and calculate K if needed (case K=0)*/
double **spk_algo(double **lnorm, int N, int *K)
{ /* Called after steps 1-2 have been made*/
    double **jacobi_output, **eigenvectors, **T, **sort_transpose;

    jacobi_output = (double **)jacobi_algo(N, lnorm);
    if (jacobi_output == NULL)
        return NULL;

    /* Transpose on eigenvectors- to make the sort easier*/
    eigenvectors = jacobi_output + 1; /* Jacobi without eigenvalues*/
    transpose(eigenvectors, N);

    /* in sort_transpose jacobi_output is being freed, there is no use in it again!*/
    sort_transpose = sort_matrix_values(jacobi_output, N);
    if (sort_transpose == NULL)
        return NULL;

    /* The Eigengap Heuristic- was told not to handle a case where k=1*/
    if (*K == 0)
        eigengap_heuristic(sort_transpose[0], N, K);

    /* Transpose on eigenvectors- get them to the right shape (vector in columns)*/
    eigenvectors = sort_transpose + 1; /* sort_transpose without eigenvalues*/
    transpose(eigenvectors, N);

    /* eigenvectors points to the start of eigenvectors, we will use only the first K vectors (first K columns) as U
     * and update eigenvectors (by renormalizing each of U’s rows) to be T */
    T = set_T(eigenvectors, N, *K);
    free_memory(sort_transpose, N + 1);

    return T;
}

/* Receives a jacobi's matrix
 * Sort first row and rows 1 to N according to the eigenvalues in first row */
double **sort_matrix_values(double **mat, int N)
{
    int i, j, max_index;
    double max_value;
    double **sort_mat = matrix_allocation(N + 1, N);
    if (sort_mat == NULL)
    {
        free_memory(mat, N + 1);
        return NULL;
    }
    for (i = 0; i < N; i++)
    {
        max_index = -1;
        max_value = -1;
        for (j = 0; j < N; j++)
        {
            /* Found new max*/
            if (max_value < mat[0][j])
            {
                max_index = j;
                max_value = mat[0][j];
            }
        }
        /* Place the i'th eigenvalue in the i'th cell and it's correspoond eigenvectors in line number i+1 */
        sort_mat[0][i] = max_value;
        sort_mat[i + 1] = mat[max_index + 1];
        mat[0][max_index] = -1;
    }
    /* free (mat=jacobi_output) */
    free(mat[0]);
    free(mat);
    return sort_mat;
}

/* Receives U (created by largest eigenvectors of jacobi), N- number of rows, K- number of columns
 * Returns T- by renormalizing each of U’s rows to have unit length */
double **set_T(double **U, int N, int K)
{
    int i, j, q;
    double sum = 0;

    double **T = matrix_allocation(N, K);
    if (T == NULL)
        return NULL;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K; j++)
        {
            /* Calculate sum once for each new row!*/
            if (j == 0)
            {
                sum = 0;
                for (q = 0; q < K; q++)
                {
                    sum += pow(U[i][q], 2);
                }
            }
            T[i][j] = (sum != 0) ? (U[i][j] / sqrt(sum)) : 0;
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
    double sum;

    double **diag_mat = matrix_allocation(N, N);
    if (diag_mat == NULL)
        return NULL;
    
    sum=0;
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
    free_memory(mul1, N);
    if (mul2 == NULL)
        return NULL;

    lnorm = I_matrix(N);
    if (lnorm == NULL)
    {
        free_memory(mul2, N);
        return NULL;
    }
    calc_sub(N, lnorm, mul2);
    free_memory(mul2, N);
    return lnorm;
}

/* Receives D- diagonal degree matrix, N- number of rows/columns
 * Updates D to D^(-0.5)*/
void cal_D12(double **diag_mat, int N)
{
    int i;
    for (i = 0; i < N; i++)
    { /* NOT TO DELETE THIS TODO!!!! TODO check if Rami answered :https://moodle.tau.ac.il/mod/forum/discuss.php?d=128591*/
        /* (assumption) diag_mat[i][i]!=0 */
        diag_mat[i][i] = 1 / sqrt(diag_mat[i][i]);
    }
}

/* Receives matrices: A,B and N- number of rows/columns
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

/* Receives matrices: A,B and N- number of rows/columns
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
/* Receives symetric matrix A- size N*N
 * Returns matrix with first row= eigenvalues, next rows are the
 * corresponding eigenvectors (each column is a vector)
 * If an error occurred returns NULL*/
double **jacobi_algo(int N, double **A)
{
    int counter, i, j, return_value;
    double c, s, offA, offA1; /*pivot element, s,c and offA and offA'*/
    /* V=V x curr_P for each update of P (update on each loop), curr_P= P calculated in each loop,
     * jacobi_result= union of eigenvalues and eigenvectors, V_to_free- in each calculate of V a new 
     * memory is allocated so V_to_free saves previoes memory pointer*/
    double **V, **curr_P, **jacobi_result,**V_to_free;
    double *eigenvalues; 
    
    /* Memory allocations- if matrix_allocation returns NULL-
     * an error occurred - free previous allocations and return NULL */
    V = I_matrix(N); /* In first iteration, V = curr_P1*/
    if (V == NULL)
        return NULL;

    curr_P = matrix_allocation(N, N);
    if (curr_P == NULL)
    {
        free_memory(V, N);
        return NULL;
    }

    eigenvalues = malloc(N * sizeof(double)); /*len of diagonal of squared matrix (NxN) is always N*/
    if (eigenvalues == NULL)
    {
        free_memory(V, N);
        free_memory(curr_P, N);
        return NULL;
    }

    counter = 0;
    return_value=1;
    offA = EPSILON_JACOBI + 1;
    offA1 = 0;

    /* Will run up to 100 iterations if it doesn't reach convergence before + in the first iteration, convergence is not relevant*/
    while ((MAX_ITER_JACOBI > counter) && ((offA - offA1 > EPSILON_JACOBI) || (counter == 0)))
    {
        counter++;
        offA = calc_off_diag(N, A);

        /* Transform the matrix A to A' */
        find_Aij(N, A, &i, &j);                   /* Finding the index of Aij - the pivot*/
        find_c_s_t(A, i, j, &c, &s);              /* Calculating c,s with the given formulas*/
        calc_curr_P(N, curr_P, i, j, c, s);       /* Calculating P matrix with the given formula*/
        calc_A1(N, A, c, s, i, j, &return_value); /* update A according to formula of A', if an error occured- return_value=0,else return_value=1 */

        V_to_free=V;
        V = calc_mul(N, V, curr_P); /* V *= curr_P */
        free_memory(V_to_free,N);
        if (return_value == 0 || V == NULL)
        { /* An error occured */
            free(eigenvalues);
            free_memory(curr_P, N);
            return NULL;
        }
        offA1 = calc_off_diag(N, A);
    }

    get_eigenvalues_from_A1(eigenvalues, N, A);            /* Getting eigenvalues from A' diagonal!*/
    jacobi_result = jacobi_eigen_merge(N, eigenvalues, V); /* Putting eigenVectors and eigenVectors together*/

    free_memory(V, N);
    free_memory(curr_P, N);
    free(eigenvalues);

    return jacobi_result;
}

/* Receives symetric matrix A- size N*N, c, s, and i, j (of pivot- Aij) and a pointer to an int return_value
 * Updates A's i and j rows and cols according to the given formula
 * If an error occurred return_value value updates to 0- else 1*/
void calc_A1(int N, double **A, double c, double s, int i, int j, int *return_value)
{
    int r;
    /* row_col_i_j[0] will represent row/col number i, row_col_i_j[1] will represent row/col number j,
       For each r: row_col_i_j[0][r]=a'[r][i]=a'[i][r] and row_col_i_j[1][r]=a'[r][j]=a'[j][r] */
    double **row_col_i_j = matrix_allocation(2, N);

    if (row_col_i_j == NULL)
    {
        *return_value = 0;
        return;
    }

    for (r = 0; r < N; r++)
    {
        /*for row/col number i */
        row_col_i_j[0][r] = (r != i && r != j) ? (c * A[r][i] - s * A[r][j]) : ((r == i) ? (c * c * A[i][i] + s * s * A[j][j] - 2 * s * c * A[i][j]) : 0);
        row_col_i_j[1][r] = (r != i && r != j) ? (c * A[r][j] + s * A[r][i]) : ((r == j) ? (s * s * A[i][i] + c * c * A[j][j] + 2 * s * c * A[i][j]) : 0);
    }
    /* update rows i,j of A */
    for (r = 0; r < N; r++)
    {
        if (r != j && r != i)
        {
            A[r][i] = row_col_i_j[0][r]; /*Updates col i (except A[i][i] and except A[j][i])*/
            A[i][r] = A[r][i];           /*Updates row i (except A[i][i] and except A[i][j])*/
            A[j][r] = row_col_i_j[1][r]; /*Updates col j (except A[j][j] and except A[j][i])*/
            A[r][j] = A[j][r];           /*Updates row j (except A[j][j] and except A[i][j])*/
        }
    }

    A[i][j] = 0;
    A[j][i] = 0;
    A[i][i] = row_col_i_j[0][i];
    A[j][j] = row_col_i_j[1][j];

    free_memory(row_col_i_j,2);
    *return_value = 1;
}

/* Receives matrix A- size N*N
 * Returns the sum of squares of all off-diagonal elements of A*/
double calc_off_diag(int N, double **A)
{
    int i, j;
    double off_A_squared = 0;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            /*Sum square of all off-diagonal elements in the Matrix*/
            off_A_squared += (i == j) ? 0 : (A[i][j] * A[i][j]);
        }
    }
    return off_A_squared;
}

/* Receives matrix mat and N- number of rows/columns
 * Updates mat to be transposed */
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

/* Receives matrix A, N- number of rows/columns and a pointer to the numbers that represent i,j
 * Updates i and j to be the locations of the largest ABSOLUTE value (the pivot) */
void find_Aij(int N, double **A, int *i_pointer, int *j_pointer)
{
    int q, l;
    double maxElem = -DBL_MAX;

    for (q = 0; q < N; ++q)
    {
        for (l = q + 1; l < N; ++l)
        {
            /* Search for the off-diagonal element with the largest absolute value */
            if (fabs(A[q][l]) > maxElem)
            {
                maxElem = fabs(A[q][l]);
                /* Changes i,j pointers*/
                *i_pointer = q;
                *j_pointer = l;
            }
        }
    }
}

/* Receives matrix A, i,j- the location of the pivot
 * Calculates c,s,t according to the given formulas */
void find_c_s_t(double **A, int i, int j, double *cPointer, double *sPointer)
{
    double theta, t;
    double signTheta = 1;

    if (A[i][j] == 0)
    {
        /* According to Rami https://moodle.tau.ac.il/mod/forum/discuss.php?d=127831 */
        *cPointer = 1;
        *sPointer = 0;
        return;
    }

    theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
    if (theta < 0)
        signTheta = -1;
    t = (signTheta) / (fabs(theta) + sqrt((theta * theta) + 1));
    *cPointer = (1) / (sqrt(1 + (t * t)));
    *sPointer = t / sqrt(1 + (t * t));
}

/* Receives matrix curr_P, N- number of rows/columns, i,j- the location of the pivot, and c,s
 * Updates P according to the given instructions */
void calc_curr_P(int N, double **curr_P, int i, int j, double c, double s)
{ /* TODO it better? (just 1 loop because you go on the diagonal)*/
    int k, l;
    /* P[i][i]=P[j][j]=c, P[i][j]=s, P[j][i]=-s 
     * on diagonal= 1, else=0*/
    for (k = 0; k < N; ++k)
    {
        for (l = 0; l < N; ++l)
        {
            if (k == l)
                curr_P[k][l] = (k == i || l == j) ? c : 1; /* 1 on the diagonal*/
            else if (k == j && l == i)
                curr_P[k][l] = -s;
            else
                curr_P[k][l] = (k == i && l == j) ? s : 0;
        }
    }
}

/* Receives matrix A1, N- number of rows/columns and a pointer to eigenvalues's array
 * Updates eigenvalues according to the diagonal of A1 */
void get_eigenvalues_from_A1(double *eigenvalues, int N, double **A1)
{
    int i;
    for (i = 0; i < N; ++i)
        eigenvalues[i] = A1[i][i];
}

/* Receives matrix eigenVectors, N- number of rows/columns and an array of the eigenvalues
 * Returns matrix with first row- eigenValues, next rows- eigenVectors */
double **jacobi_eigen_merge(int N, double *eigenValues, double **eigenVectors)
{
    double **res = NULL;
    int i;
    res = matrix_allocation(N + 1, N);
    if (res == NULL)
    {
        return NULL;
    }

    for (i = 0; i < N; i++)
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
    return count;
}

/* Receives matrices dest_mat,src_mat and number of rows and columns
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
                j--;
            }
            fgetc(ifp);
        }
    }
    rewind(ifp);
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
        if (i != num_rows - 1)
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
    if (goal == WAM || data_output == NULL)
        return data_output;

    wam_matrix = data_output;

    /* Run DDG*/
    data_output = diagonal_matrix(wam_matrix, N);
    if (data_output == NULL || goal == DDG)
    {
        free_memory(wam_matrix, N);
        return data_output;
    }

    ddg_matrix = data_output;

    /* Run LNORM*/
    data_output = laplacian_matrix(ddg_matrix, wam_matrix, N);
    free_memory(wam_matrix, N);
    free_memory(ddg_matrix, N);
    if (data_output == NULL || goal == LNORM)
        return data_output;

    lnorm_matrix = data_output;

    /* run SPK*/
    data_output = spk_algo(lnorm_matrix, N, K);
    free_memory(lnorm_matrix, N);
    return data_output;
}

/* Receives k, goal and file_name from user
 * Calculates needed information from file and call run_goal, prints result at end*/
int main(int argc, char *argv[])
{
    char *file_name;
    int N, D, K;
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
        msg_and_exit(ERROR_TYPE, 1);
    }

    print_result(data_output, N, N, goal);

    free_memory(data_input, N);
    if (goal == jacobi_g)
        free_memory(data_output, N + 1);
    else
        free_memory(data_output, N);
    fclose(ifp);
    exit(0);
}
