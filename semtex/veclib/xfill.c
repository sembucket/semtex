/*****************************************************************************
 * xfill:   x[i] = alpha.
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>


void dfill (integer n, double alpha, double* x, integer incx)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha;
}


void ifill (integer n, integer alpha, integer* x, integer incx)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha;
}


void sfill (integer n, float alpha, float* x, integer incx)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha;
}
