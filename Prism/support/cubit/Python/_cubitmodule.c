/*
 * Python cubit module
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>

#include "Python.h"
#include "cubit/cubit.h"
#include "meshobject.h"
#include "fieldobject.h"

#ifdef IRIX64
/* kludge: for some stupid reason, the complib.sgimath library makes a ref-  *
 * erence to the following functions, and I couldn't get Python to run with- *
 * out putting these in (although they're never called).                     */
int mp_setlock_  () { return 0; }
int mp_unsetlock_() { return 0; }
#endif

PyObject *build_mesh (PyObject *self, PyObject *args) 
{
  if (!PyArg_ParseTuple(args, "")) return NULL;

  return (PyObject*) mesh_alloc();
}

PyObject *build_field (PyObject *self, PyObject *args)
{
  PyObject *mesh;

  if (!PyArg_ParseTuple(args, "O", &mesh)) return NULL;

  return (PyObject*) field_alloc(MESH_HANDLE(mesh));
}

/* ------------------------------------------------------------------------- */

static PyMethodDef cubitMethods[] = {
  {"build_mesh",  build_mesh,    METH_VARARGS},
  {"build_field", build_field,   METH_VARARGS}
};

void init_cubit() {
  PyObject *module = Py_InitModule("_cubit", cubitMethods);
  cubit_init();
}
  
