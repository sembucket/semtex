/*****************************************************************************
 * xvrecp:  y[i] = 1.0 / x[i].
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif


void dvrecp (integer n, const double* x, integer incx,
                              double* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = 1.0 / x[i*incx];
}


void svrecp (integer n, const float* x, integer incx,
                              float* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = 1.0F / x[i*incx];
}
