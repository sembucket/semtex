/*****************************************************************************
 * xvexp:  y[i] = exp (x[i]).
 *
 * $Id$
 *****************************************************************************/

#include <math.h>
#include <femdef.h>


void dvexp (integer n, const double* x, integer incx,
                             double* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = exp (x[i*incx]);
}


void svexp (integer n, const float* x, integer incx,
                             float* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

#if defined(__GNUC__) || defined(__uxp__) || defined(_SX)
  for (i = 0; i < n; i++) y[i*incy] = (float) exp  (x[i*incx]);
#else
  for (i = 0; i < n; i++) y[i*incy] =         expf (x[i*incx]);
#endif
}




