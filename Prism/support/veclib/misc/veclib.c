/*
 * VECLIB -- Vector Processing Library
 *
 * This is a collection of functions for numerical and other types of opera-
 * tions on arrays.  It mimics libraries supplied for typical vector machines
 * and allows programers to concentrate on higher-level operations rather than
 * writing loops.
 *
 * ------------------------------------------------------------------------- */

#include <stdio.h>

#include "veclib/veclib.h"
#include "veclib/blas.h"
#include "veclib/lapack.h"

char   _vlib_creg [NVREG];      /* Registers for passing arguments to */
int    _vlib_ireg [NVREG];      /* call-by-reference functions.       */
float  _vlib_sreg [NVREG];
double _vlib_dreg [NVREG];

static char *revision = "$Revision$";

void veclib() 
{
  int major;
  int minor;

  if (sscanf(revision, "%*s%d.%d", &major, &minor) != 2) {
    major = 1;
    minor = 0;
  }

  printf ("veclib: version %d.%d\n", major, minor);
}

#ifdef Linux

/* The following is required under Linux so that codes linked against the *
 * f2c library don't come up with an undefined symbol error at link time. */

int MAIN__() {
  fprintf (stderr, "veclib: caught a call to MAIN__()\n");
  return 0;
}

#endif

