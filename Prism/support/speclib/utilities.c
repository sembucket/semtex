/*
 * UTILITIES
 *
 * Copyright (c) 1994 Ronald Dean Henderson
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * -------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"

/* ------------------------------------------------------------------------ *
 * ecopy(): special version of "dcopy" with negative skips                  *
 *                                                                          *
 * This does not follow the BLAS conventions on the skip.  A negative skip  *
 * goes backward in memory from the input location ["x" or "y"], not from   *
 * the end of the array ["x" + (n-1)*xskip, "y" + (n-1)*yskip].             *
 *                                                                          *
 * The name "ecopy" implies this is an edge-of-the-element copy.            *
 * ------------------------------------------------------------------------ */

int ecopy (int n, double *x, int incx, double *y, int incy)
{
  if (incx == 1 && incy == 1)
    memcpy (y, x, n*sizeof(double));
  else 
    while (n--) {
      *y = *x;
      x += incx;
      y += incy;
    }
  return 0;
}

