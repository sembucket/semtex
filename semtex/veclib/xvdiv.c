/*****************************************************************************
 * xvdiv:   z[i] = x[i] / y[i].
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>


void dvdiv (integer n,
	    const double* x, integer incx,
	    const double* y, integer incy,
	          double* z, integer incz)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = x[i*incx] / y[i*incy];
}


void svdiv (integer n,
	    const float* x, integer incx,
	    const float* y, integer incy,
	          float* z, integer incz)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = x[i*incx] / y[i*incy];
}
