/*****************************************************************************
 * xsvvtm:  z[i] = alpha - (x[i] * y[i]).
 *
 * $Id: xsvvtm.c,v 9.1 2019/05/30 06:36:14 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dsvvtm (int_t n, double alpha,
	     const double* x, int_t incx,
	     const double* y, int_t incy,
	           double* z, int_t incz)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;
  
  for (i = 0; i < n; i++) z[i*incz] = alpha - (x[i*incx] * y[i*incy]);
}


void ssvvtm (int_t n, float alpha,
	     const float* x, int_t incx,
	     const float* y, int_t incy,
	           float* z, int_t incz)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = alpha - (x[i*incx] * y[i*incy]);  
}
