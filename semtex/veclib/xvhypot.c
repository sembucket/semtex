/*****************************************************************************
 * xvhypot:  z[i] = sqrt(SQR(x[i]) + SQR(y[i])).
 *
 * $Id$
 *****************************************************************************/

#include <math.h>
#include <femdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvhypot (integer n,
	      const double* x, integer incx,
	      const double* y, integer incy,
	            double* z, integer incz)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = hypot (x[i*incx], y[i*incy]);
}


void svhypot (integer n,
	      const float* x, integer incx,
	      const float* y, integer incy,
	            float* z, integer incz)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

#if defined(__GNUC__) || defined (__uxp__) || defined(_SX)
  for (i = 0; i < n; i++) z[i*incz] = (float) hypot  (x[i*incx], y[i*incy]);
#else
  for (i = 0; i < n; i++) z[i*incz] =         hypotf (x[i*incx], y[i*incy]);
#endif
}