/*****************************************************************************
 * xcndst:  conditional assignment:  if (y[i]) z[i] = x[i].
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dcndst (integer n,
	     const double*  x, integer incx,
	     const integer* y, integer incy,
	           double*  z, integer incz)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) if (y[i*incy]) z[i*incz] = x[i*incx];
}


void icndst (integer n,
	     const integer* x, integer incx,
	     const integer* y, integer incy,
	           integer* z, integer incz)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) if (y[i*incy]) z[i*incz] = x[i*incx];
}


void scndst (integer n,
	     const float*   x, integer incx,
	     const integer* y, integer incy,
	           float*   z, integer incz)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) if (y[i*incy]) z[i*incz] = x[i*incx];
}
