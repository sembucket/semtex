/*****************************************************************************
 * Cast integer vector to floating point.
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvfloa (integer n,
	     const integer* x, integer incx,
	           double*  y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (double) x[i*incx];
}


void svfloa (integer n,
	     const integer* x, integer incx,
	           float*   y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (float) x[i*incx];
}
