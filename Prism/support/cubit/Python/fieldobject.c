/*
 * Python Field object
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include "meshobject.h"
#include "fieldobject.h" 

/* ------------------------------------------------------------------------- */

static void field_dealloc (PyObject *self) {
#ifdef DEBUG
  fprintf(stderr, "field_dealloc: %c\n", FIELD_HANDLE(self)->type);
#endif
  Field_free(FIELD_HANDLE(self));
  PyMem_DEL(self);
}

static PyObject *field_repr (PyObject *self)
{
  int   bufsiz = BUFSIZ;
  char *buf    = malloc(bufsiz);
  char *endptr = buf + bufsiz;
  char *p      = buf;

  PyObject *res;

  Field *u = FIELD_HANDLE(self);

  const int nr  = FIELD_NR(u);
  const int ns  = FIELD_NS(u);
  const int nz  = FIELD_NZ(u);

  Element *elmt;
  int i, j, k;

#ifdef DEBUG
  fprintf (stderr, "field_repr: self=%p, handle=%p\n",  self, u);
#endif

  p += sprintf(p, "Field %c:\n", FIELD_TYPE(u));
  
  for (elmt = FIELD_HEAD(u); elmt; elmt = elmt->next) {
    const int     id   = ELEMENT_ID(elmt);
    const double *uptr = FIELD_DATA(u,id);

    p += sprintf (p, "  Element %d:\n", id);
    for (k = 0; k < nz; k++) {
      for (i = 0; i < ns; i++) {
	for (j = 0; j < nr; j++) {
	  p += sprintf(p, " %#10.6g", uptr[j + ns*(i + nr*k)]);
	}
	p += sprintf (p, "\n");

	/* Check for buffer overflow on next line of data */

	if (endptr-p < 12*nr) {
	  const int buflen = p-buf;
	  bufsiz += BUFSIZ;
	  buf     = realloc(buf, bufsiz);
	  endptr  = buf + bufsiz;
	  p       = buf + buflen;
	}
      }
      if (k != nz-1) p += sprintf(p, "  --\n");
    }
  }

  *--p = '\0'; /* don't need the last newline */

  res = PyString_FromStringAndSize(buf, p-buf);
  free(buf);
  return res;
}

static PyObject *field_scal (PyObject *self, PyObject *args)
{
  double d;
  field_object *res;

  if (!PyArg_ParseTuple(args, "d", &d)) return NULL;

  res = field_alloc(FIELD_MESH(FIELD_HANDLE(self)));
  Field_copy(FIELD_HANDLE(self), FIELD_HANDLE(res));
  Field_scal(d,FIELD_HANDLE(res));

  return (PyObject*) res;
}

static PyObject *field_shift (PyObject *self, PyObject *args)
{
  double d;
  field_object *res;

  if (!PyArg_ParseTuple(args, "d", &d)) return NULL;

  res = field_alloc(FIELD_MESH(FIELD_HANDLE(self)));
  Field_copy(FIELD_HANDLE(self), FIELD_HANDLE(res));
  Field_shift(d,FIELD_HANDLE(res));

  return (PyObject*) res;
}

static PyObject *field_parse (PyObject *self, PyObject *args)
{
  char *expr;

  if (!PyArg_ParseTuple(args, "s", &expr)) return NULL;

  Field_set(FIELD_HANDLE(self), expr);
  
  Py_INCREF(self);
  return self;
}

static PyObject *field_gradient (PyObject *self, PyObject *args)
{
  int dir;
  field_object *res;
  Field *u = FIELD_HANDLE(self);
  
  if (!PyArg_ParseTuple(args, "i", &dir)) return NULL;

  if (dir < 0 || dir > 2) {
    PyErr_SetString(PyExc_ValueError, "dir must be in the range [0,2]");
    return NULL;
  }

  res = field_alloc(FIELD_MESH(u));
  Field_gradient(u, FIELD_HANDLE(res), dir);

  return (PyObject*) res;
}

static PyObject *field_copy (PyObject *self, PyObject *args) 
{
  field_object *res;
  Field *u = FIELD_HANDLE(self);

  if (!PyArg_ParseTuple(args, "")) return NULL;

  res = field_alloc(FIELD_MESH(u));
  Field_copy(u, FIELD_HANDLE(res));
  return (PyObject*) res;
}
  
static PyObject *field_FFT (PyObject *self, PyObject *args) 
{
  int dir;
  field_object *res;
  Field *u = FIELD_HANDLE(self);

  if (!PyArg_ParseTuple(args, "i", &dir)) return NULL;

  if(abs(dir) != 1) {
    PyErr_SetString (PyExc_ValueError, "dir must be 1 or -1");
    return NULL;
  }

  res = field_alloc(FIELD_MESH(u));
  Field_FFT(u, FIELD_HANDLE(res), dir);

  return (PyObject*) res;
}

static PyObject *field_settype (PyObject *self, PyObject *args) {
  char type;
  if (!PyArg_ParseTuple(args, "c", &type)) return NULL;
  return Py_BuildValue("c", FIELD_TYPE(FIELD_HANDLE(self)) = type);
}

static PyObject *field_gettype (PyObject *self, PyObject *args) {
  if (!PyArg_ParseTuple(args, "")) return NULL;
  return Py_BuildValue("c", FIELD_TYPE(FIELD_HANDLE(self)));
}

/* ------------------------------------------------------------------------- */

#define FIELD_REDUCE_FUNCTION(py_name, cubit_name) \
PyObject* py_name(PyObject *self, PyObject *args) \
{ if (!PyArg_ParseTuple(args,"")) return NULL; \
  return Py_BuildValue("d", cubit_name(FIELD_HANDLE(self))); }
  
FIELD_REDUCE_FUNCTION(field_min,      Field_min);
FIELD_REDUCE_FUNCTION(field_max,      Field_max);
FIELD_REDUCE_FUNCTION(field_amax,     Field_amax);    
FIELD_REDUCE_FUNCTION(field_nrm2,     Field_nrm2);    
FIELD_REDUCE_FUNCTION(field_integral, Field_integral);
FIELD_REDUCE_FUNCTION(field_L2,       Field_L2);      
FIELD_REDUCE_FUNCTION(field_H1,       Field_H1);      

/* ------------------------------------------------------------------------- */

static struct PyMethodDef field_methods[] = {
  {"parse",    field_parse,     METH_VARARGS},

  {"gradient", field_gradient,  METH_VARARGS},
  {"copy",     field_copy,      METH_VARARGS},
  {"FFT",      field_FFT,       METH_VARARGS},
  {"scal",     field_scal,      METH_VARARGS},
  {"shift",    field_shift,     METH_VARARGS},

  {"min",      field_min,       METH_VARARGS},
  {"max",      field_max,       METH_VARARGS},
  {"amax",     field_amax,      METH_VARARGS},
  {"integral", field_integral,  METH_VARARGS},
  {"nrm2",     field_nrm2,      METH_VARARGS},
  {"L2",       field_L2,        METH_VARARGS},
  {"H1",       field_H1,        METH_VARARGS},

  { NULL, NULL, NULL}
};

static PyObject *field_getattr (PyObject *self, char *name) {
  return Py_FindMethod(field_methods, self, name);
}

/* ------------------------------------------------------------------------- */

#define FIELD_UNARY_FUNCTION(py_name, cubit_name) \
PyObject* py_name(PyObject *op) \
{ field_object *res = field_alloc(FIELD_MESH(FIELD_HANDLE(op)));\
  cubit_name(FIELD_HANDLE(op),FIELD_HANDLE(res));\
  return (PyObject*) res; }

FIELD_UNARY_FUNCTION(field_positive, Field_copy);
FIELD_UNARY_FUNCTION(field_negative, Field_negative);
FIELD_UNARY_FUNCTION(field_abs,      Field_abs);

#define FIELD_BINARY_FUNCTION(py_name, cubit_name) \
PyObject* py_name(PyObject *op1, PyObject *op2) \
{ field_object *res = field_alloc(FIELD_MESH(FIELD_HANDLE(op1)));\
  cubit_name(FIELD_HANDLE(op1),FIELD_HANDLE(op2),FIELD_HANDLE(res));\
  return (PyObject*) res; }
      
FIELD_BINARY_FUNCTION(field_add,  Field_add);
FIELD_BINARY_FUNCTION(field_sub,  Field_sub);
FIELD_BINARY_FUNCTION(field_mult, Field_mult);
FIELD_BINARY_FUNCTION(field_div,  Field_div);

static PyNumberMethods field_as_number = {
  (binaryfunc)  field_add,            /* nb_add */
  (binaryfunc)  field_sub,            /* nb_subtract */
  (binaryfunc)  field_mult,           /* nb_multiply */
  (binaryfunc)  field_div,            /* nb_divide */
  (binaryfunc)  0,                    /* nb_remainder */
  (binaryfunc)  0,                    /* nb_divmod */
  (ternaryfunc) 0,                    /* nb_power */
  (unaryfunc)   field_negative,       /* nb_negate */
  (unaryfunc)   field_positive,       /* nb_positive */ 
  (unaryfunc)   field_abs,            /* nb_absolute */
  (inquiry)     0,                    /* nb_nonzero */
  (unaryfunc)   0,           	      /* nb_invert */
  (binaryfunc)  0,                    /* nb_lshift */
  (binaryfunc)  0,                    /* nb_rshift */
  (binaryfunc)  0,                    /* nb_and */
  (binaryfunc)  0,                    /* nb_xor */
  (binaryfunc)  0,                    /* nb_or */
  (coercion)    0,                    /* nb_coerce */
  0,                                  /* nb_int */
  0,                                  /* nb_long */
  0,                                  /* nb_float */
  0,		                      /* nb_oct */
  0,		                      /* nb_hex */
};

static PyTypeObject field_type = {
  PyObject_HEAD_INIT(0)
  0,                                  /* ob_size */
  "field",                            /* tp_name */
  sizeof(field_object),               /* tp_basizsize */
  0,                                  /* tp_itemsize */

  (destructor)  field_dealloc,        /* tp_dealloc */
  (printfunc)   0,                    /* tp_print */
  (getattrfunc) field_getattr,        /* tp_getattr */
  (setattrfunc) 0,                    /* tp_setattr */
  (cmpfunc)     0,                    /* tp_compare */
  (reprfunc)    field_repr,           /* tp_repr */

  &field_as_number,                   /* tp_as_number */
  0,                                  /* tp_as_sequence */
  0,                                  /* tp_as_mapping */
  
  (hashfunc)    0,                    /* tp_hash */
  (ternaryfunc) 0,                    /* tp_call */
  (reprfunc)    field_repr,           /* tp_str */
};

field_object *field_alloc (Mesh *mesh) 
{
  field_object *field = PyObject_NEW(field_object, &field_type);
  Field *u;

#ifdef DEBUG
  fprintf (stderr, "field_alloc\n");
#endif

  u = Field_alloc(mesh);
  FIELD_TYPE(u) = 'u';
  FIELD_HANDLE(field) = u;

  return field;
}
