/*****************************************************************************
 * xvlog:  y[i] = log(x[i]).
 *
 * $Id$
 *****************************************************************************/

#include <math.h>
#include <femdef.h>


void dvlog (integer n, const double* x, integer incx,
                         double* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = log (x[i*incx]);
}


void svlog (integer n, const float* x, integer incx,
                             float* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

#if defined(__GNUC__) || defined(__uxp__) || defined(_SX_)
  for (i = 0; i < n; i++) y[i*incy] = (float) log  (x[i*incx]);
#else
  for (i = 0; i < n; i++) y[i*incy] =         logf (x[i*incx]);
#endif
}
