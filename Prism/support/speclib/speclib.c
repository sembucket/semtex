/*
 * speclib
 *
 * $Revision$
 *
 * Author:  R. D. Henderson
 *
 * "speclib" is a library for spectral element methods.  It provides data
 * structures and routines for working with 2D functions of two space
 * dimensions, u(x,y).  These functions can be assigned values symbolically,
 * integrated and differentiated numerically, and computed implicitly by
 * solving the partial differential equation
 *
 *          \partial_x^2 u + \partial_y^2 u - \lamba^2 u = f(x,y),
 *
 * otherwise known as the Helmholtz equation.
 *
 * Other features of speclib are documented separately.
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include "speclib/speclib.h"

int speclib_init (void) {
  manager_init();         /* initialize the symbol table manager */
  BC_init();              /* initialize the BC tags */
  return 0;
}

int speclib_exit (void) {
  return 0;
}

/* ------------------------------------------------------------------------- */

int speclib_warning (const char *fmt, ...) 
{
  va_list ap;
  va_start(ap, fmt);
  fprintf (stderr, "speclib: warning: ");
  vfprintf(stderr, fmt, ap);
  va_end  (ap);
  return 0;
}

int speclib_error (const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  fprintf (stderr, "speclib: error: ");
  vfprintf(stderr, fmt, ap);
  fprintf (stderr, "\n");
  va_end  (ap);
  exit    (-1);

  return 0;
}

/* Copy the option and parameter tables to a file */

int speclib_options (FILE *fp) {
  show_options(fp);
  return 0;
}

int speclib_params (FILE *fp) {
  show_params(fp);
  return 0;
}
