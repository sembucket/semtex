/*****************************************************************************
 * xsvvtt:  z[i] = alpha  * x[i] * y[i].
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dsvvtt (integer n, double alpha,
	     const double* x, integer incx,
	     const double* y, integer incy,
	           double* z, integer incz)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;
  
  for (i = 0; i < n; i++) z[i*incz] = alpha * x[i*incx] * y[i*incy];
}


void ssvvtt (integer n, float alpha,
	     const float* x, integer incx,
	     const float* y, integer incy,
	           float* z, integer incz)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;
  
  for (i = 0; i < n; i++) z[i*incz] = alpha * x[i*incx] * y[i*incy];
}
