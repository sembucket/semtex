/*****************************************************************************
 * xvabs:  y[i] = abs(x[i]).
 *
 * $Id$
 *****************************************************************************/

#include <math.h>
#include <femdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif


void dvabs (integer n, const double* x, integer incx,
                             double* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = fabs (x[i*incx]);
}


void ivabs (integer n, const integer* x, integer incx,
                             integer* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (x[i*incx] < 0) ? -x[i*incx] : x[i*incx];
}


void svabs (integer n, const float* x, integer incx,
                             float* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

#if defined(__GNUC__) || defined(__uxp__) || defined(_SX)
  for (i = 0; i < n; i++) y[i*incy] = (float) fabs  (x[i*incx]);
#else
  for (i = 0; i < n; i++) y[i*incy] =         fabsf (x[i*incx]);
#endif
}
