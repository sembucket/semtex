/*
 * Python Mesh object
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>

#include "meshobject.h"

/* ------------------------------------------------------------------------- */

static void mesh_dealloc (PyObject *self) {
#ifdef DEBUG
  fprintf(stderr, "mesh_dealloc: %p\n", MESH_HANDLE(self));
#endif
  Mesh_free(MESH_HANDLE(self));
  PyMem_DEL(self);
}

static PyObject *mesh_info (PyObject *self, PyObject *args) 
{
  if (!PyArg_ParseTuple(args,"")) return NULL;
  Mesh_info(MESH_HANDLE(self));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *mesh_connect (PyObject *self, PyObject *args) 
{
  if (!PyArg_ParseTuple(args,"")) return NULL;
    Mesh_connect(MESH_HANDLE(self));
  Py_INCREF(Py_None);
  return Py_None;
}


PyObject *mesh_read (PyObject *self, PyObject *args)
{
  char *fname;
  Mesh *mesh;
  FILE *fp;

  if (!PyArg_ParseTuple(args, "s", &fname))
    return NULL;

  fp = fopen(fname,"r");
  if (fp == NULL) {
    PyErr_SetString(PyExc_IOError, "mesh_import failed");
    return NULL;
  }

  printf ("Reading mesh from %s\n", fname);
  Mesh_import(MESH_HANDLE(self), fp);
  fclose(fp);

  Py_INCREF(Py_None);
  return Py_None;
}


/* ------------------------------------------------------------------------- */

static struct PyMethodDef mesh_methods[] = {
  {"connect", mesh_connect, METH_VARARGS},
  {"info",    mesh_info,    METH_VARARGS},
  {"read",    mesh_read,    METH_VARARGS},

  { NULL, NULL, NULL}
};

static PyObject *mesh_getattr (PyObject *self, char *name) {
  return Py_FindMethod(mesh_methods, self, name);
}

static PyObject *mesh_repr (PyObject *self)
{
  char buf[BUFSIZ];
  Mesh *mesh = MESH_HANDLE(self);

  sprintf (buf, "Mesh: [%d %d %d] x %d", 
	   MESH_NR(mesh), MESH_NS(mesh), MESH_NZ(mesh), MESH_NELMT(mesh));

  return PyString_FromStringAndSize(buf,strlen(buf));
}

/* ------------------------------------------------------------------------- */

PyTypeObject mesh_type = {
  PyObject_HEAD_INIT(0)
  0,                                  /* ob_size */
  "mesh",                             /* tp_name */
  sizeof(mesh_object),                /* tp_basizsize */
  0,                                  /* tp_itemsize */

  (destructor)  mesh_dealloc,         /* tp_dealloc */
  (printfunc)   0,                    /* tp_print */
  (getattrfunc) mesh_getattr,         /* tp_getattr */
  (setattrfunc) 0,                    /* tp_setattr */
  (cmpfunc)     0,                    /* tp_compare */
  (reprfunc)    mesh_repr,            /* tp_repr */

  0,                                  /* tp_as_number */
  0,                                  /* tp_as_sequence */
  0,                                  /* tp_as_mapping */
  
  (hashfunc) 0,                       /* tp_hash */
  (ternaryfunc) 0,                    /* tp_call */
  (reprfunc) 0                        /* tp_str */
};

/* ------------------------------------------------------------------------- */

/* This is the only externally visible function */

mesh_object *mesh_alloc()
{
  mesh_object *mesh = PyObject_NEW(mesh_object, &mesh_type);
  mesh->handle = Mesh_alloc();
  return mesh;
}
