#ifndef FIELD_OBJECT_H
#define FIELD_OBJECT_H

/* ------------------------------------------------------------------------- *
 * field_object                                                              *
 *                                                                           *
 * This is the wrapper around cubit's Field type.  There is only one public  *
 * function, and it allows you to allocate an empty field_object.  Everything*
 * else is accessed as a PyObject method and implemented as a static C func- *
 * tion.                                                                     *
 *                                                                           *
 * For more information see the PyMethodDef table in fieldobject.c.          *
 * ------------------------------------------------------------------------- */

#include "Python.h"
#include "cubit/field.h"
#include "cubit/mesh.h"

#define FIELD_HANDLE(obj) (((field_object*)obj)->handle)

typedef struct {
  PyObject_HEAD 
  Field* handle;
} field_object;

field_object *field_alloc (Mesh *mesh);

#endif
