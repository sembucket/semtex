/*****************************************************************************
 * y = (double) x[i]
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>


void vdble (integer n, const float *x, integer incx, double *y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (double) x[i*incx];
}
