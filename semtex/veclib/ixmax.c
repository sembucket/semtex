/*****************************************************************************
 * ixmax: index of maximum value in x.
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef.h>


integer idmax (integer n, const double* x, integer incx)
{
  register integer i, imax;
  register double  xmax;

  x += (incx<0) ? (-n+1)*incx : 0;
  xmax = x[0];
  imax = 0;

  for (i = 1; i < n; i++)
   if (x[i*incx] > xmax) {
      xmax = x[i*incx];
      imax = i;
    }

  return imax;
}


integer iimax (integer n, const integer* x, integer incx)
{
  register integer i, xmax, imax;

  x += (incx<0) ? (-n+1)*incx : 0;
  xmax = x[0];
  imax = 0;

  for (i = 1; i < n; i++)
    if (x[i*incx] > xmax) {
      xmax = x[i*incx];
      imax = i;
    }

  return imax;
}


integer ismax (integer n, const float* x, integer incx)
{
  register integer i, imax;
  register float   xmax;

  x += (incx<0) ? (-n+1)*incx : 0;
  xmax = x[0];
  imax = 0;

  for (i = 1; i < n; i++)
    if (x[i*incx] > xmax) {
      xmax = x[i*incx];
      imax = i;
    }

  return imax;
}
