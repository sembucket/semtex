/*****************************************************************************
 * ixmin: index of minimum value in x.
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef.h>


integer idmin (integer n, const double* x, integer incx)
{
  register integer i, imin;
  register double  xmin;

  x += (incx<0) ? (-n+1)*incx : 0;
  xmin = x[0];
  imin = 0;

  for (i = 1; i < n; i++)
   if (x[i*incx] < xmin) {
      xmin = x[i*incx];
      imin = i;
    }

  return imin;
}


integer iimin (integer n, const integer* x, integer incx)
{
  register integer i, xmin, imin;

  x += (incx<0) ? (-n+1)*incx : 0;
  xmin = x[0];
  imin = 0;

  for (i = 1; i < n; i++)
    if (x[i*incx] < xmin) {
      xmin = x[i*incx];
      imin = i;
    }

  return imin;
}


integer ismin (integer n, const float* x, integer incx)
{
  register integer i, imin;
  register float   xmin;

  x += (incx<0) ? (-n+1)*incx : 0;
  xmin = x[0];
  imin = 0;

  for (i = 1; i < n; i++)
    if (x[i*incx] < xmin) {
      xmin = x[i*incx];
      imin = i;
    }

  return imin;
}
