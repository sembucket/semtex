/*****************************************************************************
 * xmask:  conditional assignment:  if (y[i]) z[i] = w[i]; else z[i] = x[i].
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dmask (integer n,
	    const double*  w, integer incw,
	    const double*  x, integer incx,
	    const integer* y, integer incy,
	          double*  z, integer incz)
{
  register integer i;

  w += (incw < 0) ? (-n + 1) * incw : 0;
  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;
  z += (incz < 0) ? (-n + 1) * incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = (y[i*incy]) ? w[i*incw] : x[i*incx];
}


void imask (integer n,
	    const integer* w, integer incw,
	    const integer* x, integer incx,
	    const integer* y, integer incy,
	          integer* z, integer incz)
{
  register integer i;

  w += (incw < 0) ? (-n + 1) * incw : 0;
  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;
  z += (incz < 0) ? (-n + 1) * incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = (y[i*incy]) ? w[i*incw] : x[i*incx];
}


void smask (integer n,
	    const float*   w, integer incw,
	    const float*   x, integer incx,
	    const integer* y, integer incy,
	          float*   z, integer incz)
{
  register integer i;

  w += (incw < 0) ? (-n + 1) * incw : 0;
  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;
  z += (incz < 0) ? (-n + 1) * incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = (y[i*incy]) ? w[i*incw] : x[i*incx];
}
