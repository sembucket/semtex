/*****************************************************************************
 * xsvtsp:  y[i] = alpha * x[i] + beta.
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dsvtsp (int_t n, double alpha, double beta,
	     const double* x, int_t incx,
	           double* y, int_t incy)
{
  register int_t i;
  
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = alpha * x[i*incx] + beta;
}


void svtsp (int_t n, float alpha, float beta,
	    const float* x, int_t incx,
	          float* y, int_t incy)
{
  register int_t i;
  
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = alpha * x[i*incx] + beta;
}
