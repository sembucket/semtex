/*****************************************************************************
 * xsvvttvp:  z[i] = alpha * (w[i] * x[i]) + y[i].
 *
 * $Id: xsvvttvp.c,v 9.1 2019/05/30 06:36:14 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dsvvttvp (int_t n, const double  alpha,
	       const double* w, int_t incw,
	       const double* x, int_t incx,
	       const double* y, int_t incy,
	             double* z, int_t incz)
{
  register int_t i;

  if (incw == 1 && incx == 1 && incy == 1 && incz == 1) 
   for (i = 0; i < n; i++) 
     z[i] = alpha * (w[i] * x[i]) + y[i]; 

  else {

    w += (incw<0) ? (-n+1)*incw : 0;
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    z += (incz<0) ? (-n+1)*incz : 0;

    for (i = 0; i < n; i++)
      z[i*incz] = alpha * (w[i*incw] * x[i*incx]) + y[i*incy];
  }
}


void ssvvttvp (int_t n, const float  alpha,
	       const float* w, int_t incw,
	       const float* x, int_t incx,
	       const float* y, int_t incy,
	             float* z, int_t incz)
{
  register int_t i;

  if (incw == 1 && incx == 1 && incy == 1 && incz == 1) 
   for (i = 0; i < n; i++) 
     z[i] = alpha * (w[i] * x[i]) + y[i]; 

  else {

    w += (incw<0) ? (-n+1)*incw : 0;
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    z += (incz<0) ? (-n+1)*incz : 0;

    for (i = 0; i < n; i++)
      z[i*incz] = alpha * (w[i*incw] * x[i*incx]) + y[i*incy];
  }
}
