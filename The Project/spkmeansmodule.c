//
// Created by Zohar Mosseri on 06/07/2022.
//


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