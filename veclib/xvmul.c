/*****************************************************************************
 * xvmul:  z[i] = x[i] * y[i].
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvmul (integer n,
	    const double* x, integer incx,
	    const double* y, integer incy,
	          double* z, integer incz)
{
  register integer i;
  
  if (incx == 1 && incy == 1 && incz == 1)
    for (i = 0; i < n; i++) z[i] = x[i] * y[i];
  
  else {
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    z += (incz<0) ? (-n+1)*incz : 0;

    for (i = 0; i < n; i++) z[i*incz] = x[i*incx] * y[i*incy];
  }
}


void ivmul (integer n,
	    const integer* x, integer incx,
	    const integer* y, integer incy,
	          integer* z, integer incz)
{
  register integer i;
  
  if (incx == 1 && incy == 1 && incz == 1)
    for (i = 0; i < n; i++) z[i] = x[i] * y[i];

  else {
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    z += (incz<0) ? (-n+1)*incz : 0;

    for (i = 0; i < n; i++) z[i*incz] = x[i*incx] * y[i*incy];
  }
}


void svmul (integer n,
	    const float* x, integer incx,
	    const float* y, integer incy,
	          float* z, integer incz)
{
  register integer i;

  if (incx == 1 && incy == 1 && incz == 1)
    for (i = 0; i < n; i++) z[i] = x[i] * y[i];

  else {
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    z += (incz<0) ? (-n+1)*incz : 0;

    for (i = 0; i < n; i++) z[i*incz] = x[i*incx] * y[i*incy];
  }
}
