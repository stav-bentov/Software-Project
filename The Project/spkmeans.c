
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

/*///////////////////Jacobi////////////////////////-*/
typedef struct Tuple {
    int successful;
    double* eigenvalues;
    double** eigenVectors;
}jacobiTuple;

static PyObject* fit_jacobi(PyObject *self,PyObject *args);
static jacobiTuple jacobi(int N, int max_iter, double **A, float epsilon);
static int matrix_allocation(double **mat, int size);
static int checkConvergence(int N, double **A, double **A1, float epsilon);
static double find_Aij(int N, double **A, int* iPointer, int* jPointer);
static void find_c_s_t(double **A, double aij, int i, int j, double *cPointer, double *sPointer);
static double **calc_A1(int N, double **A,double **A1, double c, double s, int i, int j);
static double **calc_curr_P(double **curr_P, int i, int j, double c, double s);
static double **calc_V(int N, double **V, double **curr_P);
static void get_eigenvalues_from_A1(double *eigenvalues, int N, double **A1);


/*///////////////////SPK////////////////////////-*/
static PyObject* fit_spk(PyObject *self,PyObject *args);
static int kMeans(int N, int K, int max_iter, float epsilon, double **Datapoints, double **Centroids, int dimension);
static int check_euclidean_norm(double **newCentroids, double **oldCentroids, int dimension, int K,float epsilon);
static int find_cluster(double **Centroids, double *Datapoint, int dimension, int K);
static void updateOldCentroid(double **newCentroids, double **oldCentroids, int dimension, int K);


/* Gets an array to be free (pointers of pointers) and their size*/
static void free_memory(double **ArrayToFree, int size){
    int i;
    for (i = 0; i < size; i++){
        free(ArrayToFree[i]);
    }
    free(ArrayToFree);
}


static jacobiTuple jacobi(int N, int max_iter, double **A, float epsilon){
    int counter, malloc_failure_check;
    int iPointer, jPointer;
    double Aij; /*pivot element*/
    double **A1; /* A' matrix */
    double **V; /*eigenVectors*/
    double **curr_P;/*P matrix - keeps changing and (V = V x curr_P)*/
    double *eigenvalues;
    double cPointer ,sPointer;

    /*allocate memory for the matrices and if there was a failure then malloc_failure_check != 0*/
    malloc_failure_check = matrix_allocation(A1, N);
    malloc_failure_check += matrix_allocation(V, N);/*todo check if size is NxN*/
    malloc_failure_check += matrix_allocation(curr_P, N);/*todo check if size is NxN*/

    eigenvalues = malloc(N * sizeof(double));/*len of diagonal of squared matrix (NxN) is always N*/
    if (eigenvalues == NULL || malloc_failure_check != 0){
        jacobiTuple structTuple = {-1, eigenvalues, V};/*FAIL*/
        return structTuple;
    }

    while ((max_iter >= counter) && (counter == 0 || !checkConvergence(N, A, A1,epsilon))){/*todo check : even if check Convergence is true in the beginning, i want to do this while loop*/
        counter++;
        /*A = A1 todo it in a function*/
        memcpy(A,A1, N*N*sizeof(double));/*todo it in a function or maybe as a loop (just copy data of A1 to A)!*/

        Aij = find_Aij(N, A, &iPointer ,&jPointer);
        find_c_s_t(A, Aij, iPointer, jPointer, &cPointer, &sPointer);
        A1 = calc_A1(N, A, A1,cPointer,sPointer,iPointer,jPointer);/*gets A, c, s, i, j*/
        curr_P = calc_curr_P(curr_P, iPointer, jPointer, cPointer, sPointer);
        V = calc_V(N, V, curr_P); /*todo - matrix multiplication*/
    }

    get_eigenvalues_from_A1(eigenvalues, N, A1); /*getting eigenvalues from A' diagonal!*/
    /*todo - remember eigenvalues must be ordered increasingly and respecting multiplicities*/
    jacobiTuple structTuple = {0, eigenvalues, V}; /* returns 1 on success, eigenvalues, eigenvectors*/
    return structTuple;
}

static double **calc_V(int N, double **V, double **curr_P) {/*todo check about complexity - maybe can do it in place or at the end just multiply them all */
    int i,j,k;
    double** c;
    matrix_allocation(c, N);

    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            for(k=0;k<N;k++){
                c[i][j]+=V[i][k]*curr_P[k][j];
            }
        }
    }
    return c;
}

static int matrix_allocation(double **mat, int size) {
    /*allocation of memory of size (nxn) and return 1 if there was a failure!*/
    int i;
    mat = malloc((sizeof(double *)) * size);
    if (mat == NULL)
        return 1;
    for (i = 0; i < size; i++) {
        mat[i] = malloc((sizeof(double)) * (size));
        if (mat[i] == NULL)
            return 1;
    }
    return 0;
}

static int checkConvergence(int N, double **A, double **A1, float epsilon){
    /* (true) iff off(A)^2 - off(A')^2 <= epsilon */
    int i,j;
    double off_A_squared, off_A1_squared;
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++) {
            if (i != j){
                off_A_squared = off_A_squared + (pow(A[i][j], 2));
                off_A1_squared = off_A1_squared + (pow(A1[i][j], 2));
            }
        }
    }
    if (off_A_squared - off_A1_squared <= epsilon)
        return 1;

    return 0;
}

static double find_Aij(int N, double **A, int* iPointer, int* jPointer) { /* finds the off-diagonal element with the largest ABSOLUTE value*/
    int q,l;
    double maxElem = -INFINITY;/*todo change - be careful*/
    for (q = 0; q < N; ++q) {
        for (l = 0; l < N; ++l) {
            if (q != l){
                if (fabs(A[q][l]) > maxElem) {
                    maxElem = A[q][l];
                    /*changes i,j pointers*/
                    *iPointer = q;
                    *jPointer = l;
                }
            }
        }
    }
    return maxElem;/*value of Aij*/
}

static void find_c_s_t(double **A, double aij, int i, int j, double *cPointer, double *sPointer) {
    double theta, t;
    double signTheta = 1;/*todo check if sign(theta) is 0 / 1 or something else*/
    theta = (A[j][j] - A[i][i]) / (2*aij);
    if (theta < 0)
        signTheta = 0;

    t = (signTheta) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    *cPointer = (1)/(sqrt(pow(t, 2) + 1));
    *sPointer = (t) * (*cPointer);/* todo check if needed * before cPointer*/
}

static double **calc_A1(int N, double **A,double **A1, double c, double s, int i, int j) {
    for (int r = 0; r < N; ++r) {
        if ((r != i) && (r != j)){
            A1[r][i] = (c * A[r][i]) - (s * A[r][j]);
            A1[r][j] = (c * A[r][j]) + (s * A[r][i]);
        }
        A1[i][i] = (pow(c,2) * A[i][i]) + (pow(s,2) * A[j][j]) - (2*(s*c*A[i][j]));
        A1[j][j] = (pow(s,2) * A[i][i]) + (pow(c,2) * A[j][j]) + (2*(s*c*A[i][j]));
    }
    A1[i][j] = ((pow(c,2) - pow(s,2)) * A[i][j]) + (s*c*(A[i][i] - A[j][j]));/*A1[i][j] = 0*/
    /* important - TODO last row which is not clear!!!!!!!*/

    return A1;
}

static double **calc_curr_P(double **curr_P, int i, int j, double c, double s) {/*todo check, maybe this function call is not necessary - complexity wise and no return val*/

    curr_P[i][i] = c;
    curr_P[i][j] = s;
    curr_P[j][i] = -s;
    curr_P[j][j] = c;
    return curr_P;/*todo - maybe not necessary this row*/
}

static void get_eigenvalues_from_A1(double *eigenvalues, int N, double **A1) {

    for (int i = 0; i < N; ++i) {
        eigenvalues[i] = A1[i][i];
    }
    /*todo check if that's the best way it can be done - no return of eigenvalues*/
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
    int l, q, m;
    l = sizeof(returnFromJacobi.eigenVectors);
    q = sizeof(returnFromJacobi.eigenvalues);
    returned_result = PyList_New(l+1);/*first row is eigenvalues then 2nd row onwards: eigenvectors*/
    current_vector = PyList_New(sizeof(returnFromJacobi.eigenvalues));/*first row : eigenvalues*/

    /*handling first row - eigenvalues*/
    /*  returned_result[0][i] = eigenvalues[i]*/
    for (i = 0; i < q; ++i) {
        PyList_SetItem(current_vector,i,Py_BuildValue("d", returnFromJacobi.eigenvalues[i]));
    }
    PyList_SetItem(returned_result,0,Py_BuildValue("O",current_vector));

    /*handling other rows - eigenvectors*/
    /*first i+1 : eigenVector of row i in eigenVectors matrix*/
    /*  returned_result[i+1][j] = eigenVectors[i][j]*/
    for(i=0; i<l; i++){
        m = sizeof(returnFromJacobi.eigenVectors[i]);
        current_vector = PyList_New(m);
        for (j=0;j<m;j++){
            PyList_SetItem(current_vector,j,Py_BuildValue("d", returnFromJacobi.eigenVectors[i][j]));
        }
        PyList_SetItem(returned_result,i+1,Py_BuildValue("O",current_vector));
    }
    free_memory(A,N);

    return returned_result;
}
/*//////////////////////////////END OF JACOBI/////////////////////////////////////*/


/*//////////////////////////////START OF SPK/////////////////////////////////////*/
static int kMeans(int N, int K, int max_iter, float epsilon, double **Datapoints, double **Centroids, int dimension){
    /*
    i, j, counter= counters for loop iterations.
    oldCentroids= saves all centroid's vectors (before change).
    */
    int i, j, counter;
    double **oldCentroids;
    int cluster;

    counter = 0;

    oldCentroids = malloc((sizeof(double *)) * K);
    if (oldCentroids == NULL)
    {
        return FAIL;
    }
    for(i=0;i<K;i++)
    {
        oldCentroids[i] = malloc((sizeof(double)) * (dimension));
        if (oldCentroids[i] == NULL)
        {
            return FAIL;
        }
        for(j=0; j<dimension; j++)
        {
            oldCentroids[i][j]=Centroids[i][j];
        }
    }

    /* Calculate k-means:
    stop iteration when number of iteration is more then man_iter
    or when all of the centroids have changed less then epsilon*/
    while (max_iter > counter)
    {
        for (i = 0; i < N; i++)
        {
            cluster = find_cluster(Centroids, Datapoints[i], dimension, K);

            /* Update for each datapoint- what number of cluster it belongs
            and for each centroids update the number of datapoints that belong to it*/
            Datapoints[i][dimension] = cluster;
            Centroids[cluster - 1][dimension] += 1;
        }

        /* Set all centroids to zero so that updated centroids can be calculated next*/
        for (i = 0; i < K; i++)
        {
            for (j = 0; j < dimension; j++)
            {
                Centroids[i][j] = 0;
            }
        }

        /* Update centroids according to the calculations*/
        for (i = 0; i < N; i++)
        {
            cluster = Datapoints[i][dimension];
            for (j = 0; j < dimension; j++)
            {
                Centroids[cluster - 1][j] += Datapoints[i][j];
            }
        }
        for (i = 0; i < K; i++)
        {
            for (j = 0; j < dimension; j++)
            {
                /*Centroids[i][j]= sum of datapoints that belong to cluster i+1
                Centroids[i][dimension]= number of datapoints that belong to cluster i*/
                Centroids[i][j] = Centroids[i][j] / Centroids[i][dimension];
            }
            Centroids[i][dimension] = 0;
        }

        /* if all centroids changed less then epsilon -done, else- countinue*/
        if (check_euclidean_norm(Centroids, oldCentroids, dimension, K,epsilon))
        {
            break;
        }

        /* make oldcentroids be the new ones for next iteration*/
        updateOldCentroid(Centroids, oldCentroids, dimension, K);

        counter++;
    }

    free_memory(oldCentroids, K);
    return SUCCESS;
}

/* Gets the new and old centroids, return 1 if all of the centroids didn't change more then epsilon,else-0*/
static int check_euclidean_norm(double **newCentroids, double **oldCentroids, int dimension, int K,float epsilon)
{
    int i, j;
    double sum;

    /*Calculate euclidean norm for each centroid*/
    for (i = 0; i < K; i++)
    {
        sum = 0;
        for (j = 0; j < dimension; j++)
        {
            sum += pow(newCentroids[i][j] - oldCentroids[i][j], 2);
        }

        /* One centroid changed more then epsilon*/
        if (sqrt(sum) >= epsilon)
            return 0;
    }
    /* Every centroids changed less then epsilon */
    return 1;
}

/* Gets the centroids and one datapoint and return datapoint's cluster*/
static int find_cluster(double **Centroids, double *Datapoint, int dimension, int K)
{
    int i, j, cluster;
    double sum, minSum;

    cluster=0; /*Default*/

    minSum = DBL_MAX;
    for (i = 0; i < K; i++)
    {
        sum = 0;
        for (j = 0; j < dimension; j++)
        {
            sum += pow((Datapoint[j] - Centroids[i][j]), 2);
        }
        if (minSum >= sum)
        {
            minSum = sum;

            /* Cluster number i+1 because it represented by index cell i*/
            cluster = i + 1;
        }
    }

    return cluster;
}


/* Gets the updated centroids and old ones- update the old centroids*/
static void updateOldCentroid(double **newCentroids, double **oldCentroids, int dimension, int K)
{
    int i, j;
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            oldCentroids[i][j] = newCentroids[i][j];
        }
    }
}

/* Gets N,K,max_iter,epsilon,dimension, Datapoints and initial Centroids lists from python and calculate their Centroids.*/
static PyObject* fit_spk(PyObject *self,PyObject *args){
    PyObject *Datapoints_PyObject;
    PyObject *Centroids_PyObject;
    PyObject *current_datapoint;
    PyObject *current_centroid;
    PyObject *current_double;
    PyObject *returned_Centroids;

    /* args= N, K, max_iter, Datapoints_array, Centroids_array, epsilon, dimension*/
    int N,K,max_iter,dimension;
    float epsilon;
    double **Datapoints;
    double **Centroids;
    int i,j;
    int return_value;

    if (!PyArg_ParseTuple(args, "iiiOOdi", &N, &K, &max_iter, &Datapoints_PyObject, &Centroids_PyObject, &epsilon,&dimension))
    {
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }

    /* Set up Datapoints and Centroids*/
    Datapoints = malloc((sizeof(double *)) * N);
    if(Datapoints==NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }
    Centroids = malloc((sizeof(double *)) * K);
    if(Centroids==NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }

    for (i = 0; i < N; i++)
    {
        Datapoints[i] = malloc((sizeof(double)) * (dimension + 1));
        if (Datapoints[i] == NULL)
        {
            PyErr_SetString(PyExc_RuntimeError, ERROR);
            return NULL;
        }
        current_datapoint = PyList_GetItem(Datapoints_PyObject, i);
        current_centroid=NULL; /*Default*/

        if(i<K)
        {
            Centroids[i] = malloc((sizeof(double)) * (dimension + 1));
            if (Centroids[i] == NULL)
            {
                PyErr_SetString(PyExc_RuntimeError, ERROR);
                return NULL;
            }
            current_centroid = PyList_GetItem(Centroids_PyObject, i);
        }

        /* Set Datapoint's and Centroid's vectors*/
        for(j=0; j<dimension; j++)
        {
            current_double=PyList_GetItem(current_datapoint,j);
            Datapoints[i][j]=PyFloat_AsDouble(current_double);

            if(i<K)
            {
                current_double=PyList_GetItem(current_centroid,j);
                Centroids[i][j]=PyFloat_AsDouble(current_double);
            }
        }

        /* Zero in last cell [dimension]*/
        Datapoints[i][j] = 0;
        if (i < K)
        {
            Centroids[i][j] = 0;
        }
    }

    /* update Centroids- is not fail then no error occured, Centroids have been updated*/
    return_value=kMeans(N,K,max_iter,epsilon,Datapoints,Centroids,dimension);
    if(return_value==FAIL)
    {
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }

    /* Convert Centroids to an array list (python)*/
    returned_Centroids=PyList_New(K);
    for(i=0; i<K; i++)
    {
        current_centroid=PyList_New(dimension);
        for (j=0;j<dimension;j++)
        {
            PyList_SetItem(current_centroid,j,Py_BuildValue("d",Centroids[i][j]));
        }
        PyList_SetItem(returned_Centroids,i,Py_BuildValue("O",current_centroid));
    }
    free_memory(Datapoints,N);
    free_memory(Centroids,K);
    return returned_Centroids;
}

/*//////////////////////////////END OF SPK/////////////////////////////////////*/


























/*static PyMethodDef Methods[]={
        {"fit_jacobi",
                (PyCFunction) fit_jacobi,
                METH_VARARGS,
                     NULL},
        {NULL,NULL,0,NULL}
};

static struct PyModuleDef moudledef={
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        NULL,
        -1,
        Methods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m=PyModule_Create(&moudledef);
    if(!m)
    {
        return NULL;
    }
    return m;
}*/


