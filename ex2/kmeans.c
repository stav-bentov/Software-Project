#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

#define INVALID "Invalid Input!"
#define ERROR "An Error Has Occurred"
#define SUCCESS 0
#define FAIL -1

/* TODO update functions- Done*/
/* TODO function description*/
static int kMeans(int N, int K, int max_iter, float epsilon, double **Datapoints, double **Centroids, int dimension);
static int check_euclidean_norm(double **newCentroids, double **oldCentroids, int dimension, int K,float epsilon);
static int find_cluster(double **Centroids, double *Datapoint, int dimension, int K);
static void free_memory(double **ArrayToFree, int size);
static void updateOldCentroid(double **newCentroids, double **oldCentroids, int dimension, int K);
static PyObject* fit(PyObject *self,PyObject *args);

static int kMeans(int N, int K, int max_iter, float epsilon, double **Datapoints, double **Centroids, int dimension)
{
    /*
    i, j, counter= counters for loop iterations.
    oldCentroids= saves all centroid's vectors (before change).
    */
    int i, j, counter;
    double **oldCentroids;
    int cluster;

    counter = 0;

    /* TODO check in python-DONE
    if (K > N || K <= 0)
    {
        printf(INVALID);
        exit(1);
    }*/

    oldCentroids = malloc((sizeof(double *)) * K);
    if (oldCentroids == NULL)
    {
        /*printf(ERROR);
        exit(1);*/
        return FAIL;
    }
    for(i=0;i<K;i++)
    {
        oldCentroids[i] = malloc((sizeof(double)) * (dimension));
        if (oldCentroids[i] == NULL)
        {
            /*printf(ERROR);
            exit(1);*/
            return FAIL;
        }
        for(j=0; j<dimension; j++)
        {
            /* TODO check if it's a variable or pointer equal*/
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

    /* TODO check if return is needed for Centroids*/

    free_memory(oldCentroids, K);
    return SUCCESS;
}

/* gets the new and old centroids, return 1 if all of the centroids didn't change more then epsilon,else-0*/
static int check_euclidean_norm(double **newCentroids, double **oldCentroids, int dimension, int K,float epsilon)
{
    int i, j;
    double sum;

    /*calculate euclidean norm for each centroid*/
    for (i = 0; i < K; i++)
    {
        sum = 0;
        for (j = 0; j < dimension; j++)
        {
            sum += pow(newCentroids[i][j] - oldCentroids[i][j], 2);
        }

        /* one centroid changed more then epsilon*/
        if (sqrt(sum) >= epsilon)
            return 0;
    }
    /* every centroids changed less then epsilon */
    return 1;
}

/* gets the centroids and one datapoint and return datapoint's cluster*/
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

            /* cluster number i+1 because it represented by index cell i*/
            cluster = i + 1;
        }
    }

    return cluster;
}

/* get all 3 arrays and free all memory that was in use*/
static void free_memory(double **ArrayToFree, int size)
{
    int i;
    for (i = 0; i < size; i++)
    {
        free(ArrayToFree[i]);
    }
    free(ArrayToFree);
}

/* gets the updated centroids and old ones- update the old centroids*/
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

static PyObject* fit(PyObject *self,PyObject *args)
{
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

    /* TODO care of ERROR */
    /* TODO check string format- \n is OK?*/
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

        /* Set datapoint's and centroid's vectors*/
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

    return_value=kMeans(N,K,max_iter,epsilon,Datapoints,Centroids,dimension);
    if(return_value==FAIL)
    {
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }
    /* TODO centroids here are good*/

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

static PyMethodDef Methods[]={
    {"fit",
    (PyCFunction) fit,
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
}