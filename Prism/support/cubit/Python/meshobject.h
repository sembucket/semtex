#ifndef MESH_OBJECT_H
#define MESH_OBJECT_H

/* ------------------------------------------------------------------------- *
 * mesh_object                                                               *
 *                                                                           *
 * This is the wrapper around cubit's Mesh type.  There is only one public   *
 * function, and it allows you to allocate an empty mesh_object.  Everything *
 * else is accessed as a PyObject method and implemented as a static C func- *
 * tion.                                                                     *
 *                                                                           *
 * For more information see the PyMethodDef table in meshobject.c.           *
 * ------------------------------------------------------------------------- */

#include "Python.h"
#include "cubit/mesh.h"

#define MESH_HANDLE(obj) (((mesh_object*)obj)->handle)

typedef struct {
  PyObject_HEAD 
  Mesh* handle;
} mesh_object;

mesh_object *mesh_alloc();

#endif
