/*****************************************************************************
 * xvneg:  y[i] = -x[i].
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>


void dvneg (integer n, const double* x, integer incx,
                             double* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = -x[i*incx];
}


void ivneg (integer n, const integer* x, integer incx,
                             integer* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = -x[i*incx];
}


void svneg (integer n, const float* x, integer incx,
                             float* y, integer incy)
{
  register integer i;
  
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = -x[i*incx];
}
