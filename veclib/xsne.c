/*****************************************************************************
 * xsne:  y[i] = alpha != x[i].
 *
 * $Id: xsne.c,v 9.1 2019/05/30 06:36:14 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dsne (int_t n, double alpha,
	   const double* x, int_t incx,
	         double* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (alpha != x[i*incx]) ? 1 : 0;
}


void isne (int_t n, int_t alpha,
	   const int_t* x, int_t incx,
	         int_t* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (alpha != x[i*incx]) ? 1 : 0;
}


void ssne (int_t n, float alpha,
	   const float* x, int_t incx,
	         float* y, int_t incy)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (alpha != x[i*incx]) ? 1 : 0;
}
