/*****************************************************************************
 * xsle:  y[i] = alpha <= x[i].
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dsle (integer n, double alpha,
	   const double* x, integer incx,
	         double* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (alpha <= x[i*incx]) ? 1 : 0;
}


void isle (integer n, integer alpha,
	   const integer* x, integer incx,
	         integer* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (alpha <= x[i*incx]) ? 1 : 0;
}


void ssle (integer n, float alpha,
	   const float* x, integer incx,
	         float* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (alpha <= x[i*incx]) ? 1 : 0;
}
