/*****************************************************************************
 * xclip: clip to interval [alpha,beta]: if y[i] = MIN(MAX(x[i],alpha),beta).
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

#define MIN(a, b)  ((a) < (b) ?     (a) : (b))
#define MAX(a, b)  ((a) > (b) ?     (a) : (b))

void dclip (integer n, const double alpha, const double beta,
	    const double* x, integer incx,
	          double* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MIN(MAX(x[i*incx],alpha),beta);
}


void iclip (integer n, const integer alpha, const integer beta,
	    const integer* x, integer incx,
	          integer* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MIN(MAX(x[i*incx],alpha),beta);
}


void sclip (integer n, const float alpha, const float beta,
	    const float* x, integer incx,
	          float* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MIN(MAX(x[i*incx],alpha),beta);
}
