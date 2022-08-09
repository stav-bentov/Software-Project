#include <Python.h>
#include "spkmeans.h"

static PyObject *fit(PyObject *self, PyObject *args)
{
    PyObject *Datapoints_PyObject;
    PyObject *Centroids_PyObject;
    PyObject *current_datapoint;
    PyObject *current_double;
    PyObject *returned_result;
    PyObject *current_vector;
    PyObject *current_centroid;

    /* (Goal= wam, ddg, lnorm, jacobi, spk(1)): args= N, K, D, Datapoints/matrix, goal */
    /* (Goal= spk(2)): args= N, K, D, Datapoints/matrix, goal, Centroids*/
    /* args= N, K, D, Datapoints/matrix, goal, Centroids*/
    int N, K, D, i, j, rows, cols, cols_allocation, return_value;
    double **Datapoints, **Centroids;
    enum Goal goal;
    double **goal_result;

    /*receiving args from Python program*/
    if (!PyArg_ParseTuple(args, "iiiOiO", &N, &K, &D, &Datapoints_PyObject, &goal, &Centroids_PyObject))
    { /*todo check if goal sent from python is a number*/
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }

    /* only in spk_ex2- when we return to fit in the second time, Datapoints need to have D+1 cols*/
    if(goal==SPK_EX2)
        cols_allocation=D+1;
    else
        cols_allocation=D;

    /* Set up Datapoints and Centroids's matrix*/
    Datapoints = matrix_allocation(N, cols_allocation);
    if (Datapoints == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }
    current_centroid=NULL;
    Centroids=NULL;
    /* Set Centroids only in case spk_ex2!*/
    if(goal == SPK_EX2)
    {
        Centroids = matrix_allocation(K, cols_allocation);
        if (Centroids == NULL)
        {
            free_memory(Datapoints, N);
            PyErr_SetString(PyExc_RuntimeError, ERROR);
            return NULL;
        }
    }
    
    /* Fill matrix value as given matrix*/
    for (i = 0; i < N; i++)
    {
        current_datapoint = PyList_GetItem(Datapoints_PyObject, i);
        if (i < K && goal == SPK_EX2)
            current_centroid = PyList_GetItem(Centroids_PyObject, i);

        /*Set up each of vector*/
        for (j = 0; j < D; j++)
        {
            current_double = PyList_GetItem(current_datapoint, j);
            Datapoints[i][j] = PyFloat_AsDouble(current_double);
            if (i < K && goal == SPK_EX2)
            {
                current_double = PyList_GetItem(current_centroid, j);
                Centroids[i][j] = PyFloat_AsDouble(current_double);
            }
        }

        /* Only in spk_ex2- Zero in last cell [dimension]*/
        if(goal == SPK_EX2)
        {
            Datapoints[i][j] = 0;
            if (i < K)
            {
                Centroids[i][j] = 0;
            }
        }
    }

    if (goal == SPK_EX2)
    {
        return_value = kMeans(N, K, Datapoints, Centroids, D);
        if (return_value == FAIL)
        { /*todo: add free_mem*/
            PyErr_SetString(PyExc_RuntimeError, ERROR);
            return NULL;
        }
        goal_result = Centroids;
        rows = K; /* cols=number of centroids=K*/
        cols=D; /*cols= dimension*/
    }
    else
    {
        goal_result = run_goal(goal, Datapoints, N, D, &K);
        if (goal_result == NULL) /*todo check if that's how Stav wants it to be*/
        {
            free_memory(Datapoints, N);
            PyErr_SetString(PyExc_RuntimeError, ERROR);
            return NULL;
        }
        rows = (goal == JACOBI) ? (N + 1) : N; /*jacobi needs N+1 rows, rows= Number of points (int this case)*/
        cols=N; /* in wam,ddg,lnorm,jacobi*/
        if(goal==SPK)
            cols=K; /* in spk- T dimesnion are n*k (updated or original k) */
    }

    /* Convert result_matrix to an array list (python)*/
    returned_result = PyList_New(rows);
    for (i = 0; i < rows; ++i)
    {
        current_vector = PyList_New(cols);
        for (j = 0; j < cols; j++)
        {
            PyList_SetItem(current_vector, j, Py_BuildValue("d", goal_result[i][j]));
        }
        PyList_SetItem(returned_result, i, Py_BuildValue("O", current_vector));
    }
    free_memory(Datapoints, N);
    free_memory(goal_result, rows);

    return returned_result;
}

static PyMethodDef Methods[] = {
        {"fit",
                (PyCFunction)fit,
                     METH_VARARGS,
                        NULL},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moudledef = {
        PyModuleDef_HEAD_INIT,
        "my_spkmeans",
        NULL,
        -1,
        Methods
};

PyMODINIT_FUNC PyInit_my_spkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&moudledef);
    if (!m)
    {
        return NULL;
    }
    return m;
}
