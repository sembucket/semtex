/*****************************************************************************
 * xvvtvvtm:  z[i] = (v[i] * w[i]) - (x[i] * y[i]).
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvvtvvtm (integer n,
	       const double* v, integer incv,
	       const double* w, integer incw,
	       const double* x, integer incx,
	       const double* y, integer incy,
	             double* z, integer incz)
{
  register integer i;

  if (incv == 1 && incw == 1 && incx == 1 && incy == 1 && incz == 1) 
   for (i = 0; i < n; i++) 
     z[i] = v[i] * w[i] - x[i] * y[i]; 

  else {

    w += (incw<0) ? (-n+1)*incw : 0;
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    z += (incz<0) ? (-n+1)*incz : 0;

    for (i = 0; i < n; i++)
      z[i*incz] = v[i*incv] * w[i*incw] - x[i*incx] * y[i*incy];
  }
}


void svvtvvtm (integer n,
	       const float* v, integer incv,
	       const float* w, integer incw,
	       const float* x, integer incx,
	       const float* y, integer incy,
	             float* z, integer incz)
{
  register integer i;

  if (incv == 1 && incw == 1 && incx == 1 && incy == 1 && incz == 1) 
   for (i = 0; i < n; i++) 
     z[i] = v[i] * w[i] - x[i] * y[i]; 

  else {

    w += (incw<0) ? (-n+1)*incw : 0;
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    z += (incz<0) ? (-n+1)*incz : 0;

    for (i = 0; i < n; i++)
      z[i*incz] = v[i*incv] * w[i*incw] - x[i*incx] * y[i*incy];
  }
}
