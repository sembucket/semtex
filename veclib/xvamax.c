/*****************************************************************************
 * xvamax:  z[i] = MAX(ABS(x[i]), ABS(y[i])).
 *
 * $Id$
 *****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include <femdef.h>

#define MAX(x, y) ( ((x)>(y)) ? (x) : (y))


void dvamax (integer n, 
	     const double* x, integer incx,
	     const double* y, integer incy,
	           double* z, integer incz)
{
  register integer i;
  register double  absx, absy;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) {
    absx      = fabs (x[i*incx]);
    absy      = fabs (y[i*incy]);
    z[i*incz] = MAX (absx, absy);
  }
}


void ivamax (integer n,
	     const integer* x, integer incx,
	     const integer* y, integer incy,
	           integer* z, integer incz)
{
  register integer i, absx, absy;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) {
    absx      = abs (x[i*incx]);
    absy      = abs (y[i*incy]);
    z[i*incz] = MAX (absx, absy);
  }
}


void svamax (integer n,
	     const float* x, integer incx,
	     const float* y, integer incy,
	           float* z, integer incz)
{
  register integer i;
  register float   absx, absy;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) {
#if defined(__GNUC__) || defined(__uxp__) || defined(_SX)
    absx = (float) fabs (x[i*incx]);
    absy = (float) fabs (y[i*incy]);
#else
    absx = fabsf (x[i*incx]);
    absy = fabsf (y[i*incy]);
#endif
    z[i*incz] = MAX (absx, absy);
  }
}
