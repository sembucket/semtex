/*
 * cubit
 *
 * $Revision$
 *
 * Author:  R. D. Henderson
 *
 * "cubit" is a library for spectral element methods.  It provides data
 * structures and routines for working with 2D functions of two space
 * dimensions, u(x,y).  These functions can be assigned values symbolically,
 * integrated and differentiated numerically, and computed implicitly by
 * solving the partial differential equation
 *
 *          \partial_x^2 u + \partial_y^2 u + \lamba^2 u = f(x,y),
 *
 * otherwise known as the Helmholtz equation.
 *
 * "cubit" also supported adaptive mesh refinement via regular subdivisions
 * of an existing mesh.  The resulting discretization of space follows the
 * structure of a quadtree, with successive levels of refinement nested 
 * inside coarser levels of the mesh.
 *
 * ------------------------------------------------------------------------- */

#include "cubit.h"

int cubit_init() {
  manager_init();         /* initialize the symbol table manager */
  return 0;
}

int cubit_exit() {
  return 0;
}
