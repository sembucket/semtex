/*****************************************************************************
 * xclip: clip to interval [alpha,beta]: y[i] = MIN(MAX(x[i],alpha),beta).
 *
 * xclipup: clip on lower limit alpha:   y[i] = MAX(x[i],alpha).
 * xclipdn: clip on upper limit alpha:   y[i] = MIN(x[i],alpha).
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


void dclipup (integer n, const double alpha,
	      const double* x, integer incx,
	            double* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MAX(x[i*incx], alpha);
}


void iclipup (integer n, const integer alpha,
	      const integer* x, integer incx,
	            integer* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MAX(x[i*incx], alpha);
}


void sclipup (integer n, const float alpha,
	      const float* x, integer incx,
	            float* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MAX(x[i*incx], alpha);
}


void dclipdn (integer n, const double alpha,
	      const double* x, integer incx,
	            double* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MIN(x[i*incx], alpha);
}


void iclipdn (integer n, const integer alpha,
	      const integer* x, integer incx,
	            integer* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MIN(x[i*incx], alpha);
}


void sclipdn (integer n, const float alpha,
	      const float* x, integer incx,
	            float* y, integer incy)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = MIN(x[i*incx], alpha);
}
