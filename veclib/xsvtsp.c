/*****************************************************************************
 * xsvtsp:  y[i] = alpha * x[i] + beta.
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dsvtsp (integer n, double alpha, double beta,
	     const double* x, integer incx,
	           double* y, integer incy)
{
  register integer i;
  
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = alpha * x[i*incx] + beta;
}


void svtsp (integer n, float alpha, float beta,
	    const float* x, integer incx,
	          float* y, integer incy)
{
  register integer i;
  
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = alpha * x[i*incx] + beta;
}
