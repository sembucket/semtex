/*****************************************************************************
 * xvsin:  y[i] = sin(x[i]).
 *
 * $Id$
 *****************************************************************************/

#include <math.h>
#include <cfemdef>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvsin (integer n, const double* x, integer incx,
                             double* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = sin (x[i*incx]);
}


void svsin (integer n, const float* x, integer incx,
                             float* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

#if defined(__GNUC__) || defined(__uxp__) || defined(_SX)
  for (i = 0; i < n; i++) y[i*incy] = (float) sin  (x[i*incx]);
#else
  for (i = 0; i < n; i++) y[i*incy] =         sinf (x[i*incx]);
#endif
}
