/*****************************************************************************
 * y[i] = (float) x[i]
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef>



void vsngl (integer n, const double *x, integer incx, float *y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (float) x[i*incx];
}
