
static PyObject* fit(PyObject *self,PyObject *args){

    PyObject *Datapoints_PyObject; /*A matrix*/
    PyObject *current_datapoint;
    PyObject *current_double;
    PyObject *returned_result;
    PyObject *current_vector;

    /* args= N, K, max_iter, A, epsilon*/
    int N,K,D,i,j,rows;
    double **Datapoints;
    enum Goal goal;
    double ** goal_result;
    /*receiving args from Python program*/
    if (!PyArg_ParseTuple(args, "iiiOi", &N, &K, &D, &Datapoints_PyObject, &goal)){/*todo check if goal sent from python is a number*/
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }
    /* Set up Datapoints matrix*/
    Datapoints = matrix_allocation(N, D); /* todo check if cols == D or N*/
    if(Datapoints == NULL){
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }
    for (i = 0; i < N; i++){
        current_datapoint = PyList_GetItem(Datapoints_PyObject, i);
        /*Set up each of Datapoints vectors*/
        for(j = 0; j < D; j++){
            current_double=PyList_GetItem(current_datapoint,j);
            Datapoints[i][j]=PyFloat_AsDouble(current_double);
        }
    }
    rows = (goal == 4) ? N+1 : N;/*jacobi needs N+1 rows*/
    goal_result = matrix_allocation(rows, D);
    goal_result = run_goal(goal, Datapoints, N, D, K);
    if(goal_result == NULL) /*todo check if thats how Stav wants it to be*/
    {
        PyErr_SetString(PyExc_RuntimeError, ERROR);
        return NULL;
    }

    /* Convert result_matrix to an array list (python)*/

    returned_result = PyList_New(rows);
    for (i = 0; i < rows; ++i)
    {
        current_vector = PyList_New(D);
        for (j = 0; j < D; j++)
        {
            PyList_SetItem(current_vector, j, Py_BuildValue("d", goal_result[i][j]));
        }
        PyList_SetItem(returned_result,i,Py_BuildValue("O",current_vector));
    }
    free_memory(Datapoints, N);
    free_memory(goal_result, rows);
    return returned_result;
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