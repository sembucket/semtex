/*****************************************************************************
 * xsvvttvp:  z[i] = alpha * (w[i] * x[i]) + y[i].
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dsvvttvp (integer n, const double  alpha,
	       const double* w, integer incw,
	       const double* x, integer incx,
	       const double* y, integer incy,
	             double* z, integer incz)
{
  register integer i;

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


void ssvvttvp (integer n, const float  alpha,
	       const float* w, integer incw,
	       const float* x, integer incx,
	       const float* y, integer incy,
	             float* z, integer incz)
{
  register integer i;

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
