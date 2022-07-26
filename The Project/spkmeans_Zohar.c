
#define PY_SSIZE_T_CLEAN
#include <Python/Python.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <float.h>


#define INVALID "Invalid Input!"
#define ERROR "An Error Has Occurred"
#define SUCCESS 0
#define FAIL -1


/*todo! take all fit functions to h file and to module file!!!!!*/

static void free_memory(double **ArrayToFree, int size);

/* Gets an array to be free (pointers of pointers) and their size*/
static void free_memory(double **ArrayToFree, int size){
    int i;
    for (i = 0; i < size; i++){
        free(ArrayToFree[i]);
    }
    free(ArrayToFree);
}

/*///////////////////Jacobi////////////////////////-*/
typedef struct Tuple {
    int successful;
    double* eigenvalues;
    double** eigenVectors;
}jacobiTuple;

static PyObject* fit_jacobi(PyObject *self,PyObject *args);
static jacobiTuple jacobi(int N, int max_iter, double **A, float epsilon);
void identity_matrix(double **mat, int N);
static double ** matrix_allocation(double **mat, int rows, int columns);
static void matrixcopy(int N, double ** destmat, double ** srcmat);
static int checkConvergence(int N, double **A, double **A1, float epsilon);
static double find_Aij(int N, double **A, int* iPointer, int* jPointer);
static void find_c_s_t(double **A, double aij, int i, int j, double *cPointer, double *sPointer);
/*static void calc_A1(int N, double **A,double **A1, double c, double s, int i, int j);*/
static void calc_curr_P(int N, double **curr_P, int i, int j, double c, double s);
static void matrix_multiplication(int N, double **src1, double **src2, double **dst);
static void get_eigenvalues_from_A1(double *eigenvalues, int N, double **A1);
static void printmatrices(int rows, int columns, double **mat_to_print);
static void printarray(int N, double *arr);
static void transpose(double **mat, int N);




static jacobiTuple jacobi(int N, int max_iter, double **A, float epsilon){
    int counter = 0;
    int iPointer, jPointer;
    double Aij; /*pivot element*/
    double **A1 = NULL; /* A' matrix */
    double **V = NULL, **V_tag=  NULL; /*eigenVectors*/
    double **curr_P =  NULL;/*P matrix - keeps changing and (V = V x curr_P)*/
    double *eigenvalues;
    double cPointer ,sPointer;

    A1 = matrix_allocation(A1, N, N);
    V = matrix_allocation(V, N, N);
    identity_matrix(V, N); /*V equals identity matrix at the beginning*/
    curr_P = matrix_allocation(curr_P, N, N);
    eigenvalues = malloc(N * sizeof(double));/*len of diagonal of squared matrix (NxN) is always N*/
    if ((eigenvalues == NULL) || (A1 == NULL) || (V == NULL) || (curr_P == NULL) ){
        jacobiTuple structTuple = {FAIL, eigenvalues, V};/*FAIL*/
        return structTuple;
    }


    /*todo check*/
    while ((max_iter > counter) && (!checkConvergence(N, A, A1,epsilon) || (counter == 0) )){/*todo check : even if check Convergence is true in the beginning, i want to do this while loop*/
        counter++;
        /*A = A1*/
        if (counter != 1)
            matrixcopy(N, A, A1);

        Aij = find_Aij(N, A, &iPointer ,&jPointer);
        find_c_s_t(A, Aij, iPointer, jPointer, &cPointer, &sPointer);
        calc_curr_P(N, curr_P, iPointer, jPointer, cPointer, sPointer);
        /*calc_A1(N, A, A1,cPointer,sPointer,iPointer,jPointer);gets A, c, s, i, j*/
        transpose(curr_P, N);/* P -> P_transpose */
        matrix_multiplication(N, curr_P, A, A1); /* A' = P_transpose*A */
        transpose(curr_P, N);/* P_transpose -> P */
        matrix_multiplication(N, A1, curr_P, A1); /* A' = (P_transpose*A)*P */
        matrix_multiplication(N, V, curr_P, V);/* V = V * curr_P*/
    }

    get_eigenvalues_from_A1(eigenvalues, N, A1); /*getting eigenvalues from A' diagonal!*/
    printarray(N,  eigenvalues);
    /*printmatrices(N, N, V);*/
    /*todo - remember eigenvalues must be ordered increasingly and respecting multiplicities*/
    jacobiTuple structTuple = {SUCCESS, eigenvalues, V}; /* returns 1 on success, eigenvalues, eigenvectors*/
    return structTuple;
}

static void transpose(double **mat, int N) {
    int i,j;
    double tmp;

    for(i=0; i<N; i++){
        for(j=i+1; j<N; j++){
            tmp = mat[i][j];
            mat[i][j] = mat[j][i];
            mat[j][i] = tmp;
        }
    }
}


void identity_matrix(double **mat, int N){
    int i,j;

    for(i=0;i<N;i++) {
        for (j = 0; j < N; j++) {
            mat[i][j] = i == j ? 1 : 0;
        }
    }
}

static void printarray(int N, double *arr) {
    for (int i = 0; i < N; ++i) {
        printf("arr[%d] = %lf \n", i,arr[i]);
    }
}

static void printmatrices(int rows, int columns, double **mat_to_print) {
    int i, j;

    for (i=0; i<rows; i++) {
        for (j = 0; j < columns; j++) {
            printf("mat[%d][%d] = %lf ", i,j,mat_to_print[i][j]);
        }
        printf("\n");
    }
}


static void matrixcopy(int N, double ** destmat, double ** srcmat){
    int i, j;

    for (i=0; i<N; i++) { /* rad-nr */
        for (j = 0; j < N; j++) { /* kolumn-nr */
            destmat[i][j] = srcmat[i][j];
        }
    }
}


static void matrix_multiplication(int N, double **src1, double **src2, double **dst){
    double **temp;
    int i,j,k;
    temp = matrix_allocation(temp, N, N);

    for(i=0; i<N; i++){ /* temp = src1*src2 */
        for(j=0; j<N; j++){
            temp[i][j] = 0;
            for(k=0; k<N; k++)
                temp[i][j]+= src1[i][k] * src2[k][j];
        }
    }
    for(i=0; i<N; i++) { /* dst = temp */
        for (j = 0; j < N; j++) {
            dst[i][j] = temp[i][j];
        }
    }

    free_memory(temp, N);
}

static double ** matrix_allocation(double **mat, int rows, int columns) { /*good*/
    /*allocation of memory of size (nxn) and return 1 if there was a failure!*/
    int i;
    mat = malloc((sizeof(double *)) * rows);
    if (mat == NULL)
        return NULL;
    for (i = 0; i < rows; i++) {
        mat[i] = malloc((sizeof(double)) * (columns));
        if (mat[i] == NULL)
            return NULL;
    }
    return mat;
}

static int checkConvergence(int N, double **A, double **A1, float epsilon){
    /* (true) iff off(A)^2 - off(A')^2 <= epsilon */
    int i,j;
    double off_A_squared = 0, off_A1_squared = 0;

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++) {
            off_A_squared += i==j ? 0 : pow(A[i][j],2);
            off_A1_squared += i==j ? 0 : pow(A1[i][j],2);
        }
    }
    printf("covr : %lf\n", off_A_squared - off_A1_squared);
    if (off_A_squared - off_A1_squared <= epsilon)
        return 1;
    return 0;
}

static double find_Aij(int N, double **A, int* iPointer, int* jPointer) { /* finds the off-diagonal element with the largest ABSOLUTE value*/
    int q,l;
    double maxElem = -DBL_MAX;/*todo change - be careful*/

    for (q = 0; q < N; ++q) {
        for (l = 0; l < N; ++l) {
            if ((q != l) && fabs(A[q][l]) > maxElem) {
                maxElem = fabs(A[q][l]);
                /*changes i,j pointers*/
                *iPointer = q;
                *jPointer = l;
            }
        }
    }
    return maxElem;/*value of Aij*/
}

static void find_c_s_t(double **A, double aij, int i, int j, double *cPointer, double *sPointer) {
    double theta, t;
    double signTheta = 1;/*todo check if sign(theta) is 0 / 1 or something else*/
    theta = (A[j][j] - A[i][i]) / (2*A[i][j]);

    if (theta < 0)
        signTheta = -1;
    t = (signTheta) / (fabs(theta) + sqrt(theta*theta + 1));
    *cPointer = (1)/(sqrt(1 + t*t));
    *sPointer = t/sqrt(1+t*t);
}

/*static void calc_A1(int N, double **A,double **A1, double c, double s, int i, int j) {
    int r;
    for (r = 0; r < N; ++r) {
        if ((r != i) && (r != j)){
            A1[r][i] = (c * A[r][i]) - (s * A[r][j]);
            A1[r][j] = (c * A[r][j]) + (s * A[r][i]);
        }
    }
    A1[i][i] = (pow(c,2) * A[i][i]) + (pow(s,2) * A[j][j]) - (2*(s*c*A[i][j]));
    A1[j][j] = (pow(s,2) * A[i][i]) + (pow(c,2) * A[j][j]) + (2*(s*c*A[i][j]));
    A1[i][j] = ((pow(c,2) - pow(s,2)) * A[i][j]) + (s*c*(A[i][i] - A[j][j]));
    //A1[i][j] = 0;//
    // important - TODO last row which is not clear!!!!!!!
}*/

static void calc_curr_P(int N, double **curr_P, int i, int j, double c, double s) {/*todo check, maybe this function call is not necessary - complexity wise and no return val*/
    int k,l;

    for (k = 0; k < N; ++k) {
        for (l = 0; l < N; ++l) {
            if (k == l)
                curr_P[k][l] = (k == i || l == j) ? c : 1;
            else if (k == j && l == i)
                curr_P[k][l] = -s;
            else curr_P[k][l] = (k == i && l == j) ? s : 0;
        }
    }

}

static void get_eigenvalues_from_A1(double *eigenvalues, int N, double **A1) {
    int i;
    for (i = 0; i < N; ++i)
        eigenvalues[i] = A1[i][i];
}

static double** jacobi_eigen_merge(int N, double* eigenValues, double ** eigenVectors){
    double** res = NULL;
    int i,j;
    res = matrix_allocation(res, N+1 , N);

    for (i = 0; i < N; ++i)
        res[0][i] = eigenValues[i];

    for(i=0; i<N; i++) {
        for (j = 0; j < N; j++) {
            res[i+1][j] = eigenVectors[i][j];
        }
    }
    return res;
}

/* Gets N,K,max_iter,A,epsilon from python and calculate their Centroids.*/
static PyObject* fit_jacobi(PyObject *self,PyObject *args){

    PyObject *A_PyObject; /*A matrix*/
    PyObject *current_datapoint;
    PyObject *current_double;
    PyObject *returned_result;
    PyObject *current_vector;

    /* args= N, K, max_iter, A, epsilon*/
    int N,K,max_iter;
    float epsilon;
    double **A;
    int i,j;
    jacobiTuple returnFromJacobi;
    /*receiving args from Python program*/
    if (!PyArg_ParseTuple(args, "iiiOd", &N, &K, &max_iter, &A_PyObject, &epsilon)){
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }
    /* Set up A matrix*/
    A = malloc((sizeof(double *)) * N);
    if(A==NULL){
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }
    for (i = 0; i < N; i++){
        A[i] = malloc((sizeof(double)) * (N));
        if (A[i] == NULL)
        {
            PyErr_SetString(PyExc_RuntimeError, ERROR);
            return NULL;
        }
        current_datapoint = PyList_GetItem(A_PyObject, i);

        /* Set up each of A vectors*/
        for(j=0; j<N; j++){
            current_double=PyList_GetItem(current_datapoint,j);
            A[i][j]=PyFloat_AsDouble(current_double);
        }
    }
    /* returns (SUCCESS/FAIL status, eigenvalues, eigenVectors)*/
    returnFromJacobi = jacobi(N, max_iter, A, epsilon);

    if(returnFromJacobi.successful == FAIL)
    {
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }

    /* Convert Centroids to an array list (python)*/
    /*int l, q, m;
    l = sizeof(returnFromJacobi.eigenVectors);
    q = sizeof(returnFromJacobi.eigenvalues);*/
    returned_result = PyList_New(N+1);/*first row is eigenvalues then 2nd row onwards: eigenvectors*/
    current_vector = PyList_New(N);/*first row : eigenvalues*/

    /*handling first row - eigenvalues*/
    /*  returned_result[0][i] = eigenvalues[i]*/
    for (i = 0; i < N; ++i) {
        PyList_SetItem(current_vector,i,Py_BuildValue("d", returnFromJacobi.eigenvalues[i]));
    }
    PyList_SetItem(returned_result,0,Py_BuildValue("O",current_vector));

    /*handling other rows - eigenvectors*/
    /*first i+1 : eigenVector of row i in eigenVectors matrix*/
    /*  returned_result[i+1][j] = eigenVectors[i][j]*/
    for(i=0; i<N; i++){
        /*m = sizeof(returnFromJacobi.eigenVectors[i]);*/
        current_vector = PyList_New(N);
        for (j=0;j<N;j++){
            PyList_SetItem(current_vector,j,Py_BuildValue("d", returnFromJacobi.eigenVectors[i][j]));
        }
        PyList_SetItem(returned_result,i+1,Py_BuildValue("O",current_vector));
    }
    free_memory(A,N);

    return returned_result;
}
/*//////////////////////////////END OF JACOBI/////////////////////////////////////*/

int main(){
    int size = 3;
    double ** m = NULL;
    double ** result = NULL;
    jacobiTuple returnFromJacobi;
    m = matrix_allocation(m,size, size);
    m[0][0] = 0.88612627;
    m[0][1] = 0.36226522;
    m[0][2] = 0.01146638;
    m[1][0] = 0.36226522;
    m[1][1] = 0.60832993;
    m[1][2] = 0.17301408;
    m[2][0] = 0.01146638;
    m[2][1] = 0.17301408;
    m[2][2] = 0.75444075;

    returnFromJacobi = jacobi(size, 100, m, 0.00001);
    result = jacobi_eigen_merge(size, returnFromJacobi.eigenvalues, returnFromJacobi.eigenVectors);
    printmatrices(size+1, size, result);

    /*int l,a,q,j,i, p;
    l = sizeof(returnFromJacobi.eigenVectors);
    q = sizeof(returnFromJacobi.eigenvalues);
    for (p = 0; p < q; p++) {
        printf("eigenvalue: %lf \n", returnFromJacobi.eigenvalues[p]);
    }

    for(i=0; i<l; i++){
        a = sizeof(returnFromJacobi.eigenVectors[i]);
        for (j=0;j<a;j++){
            printf("eigenvector: %lf \n", returnFromJacobi.eigenVectors[i][j]);
        }
    }*/
}
