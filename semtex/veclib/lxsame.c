/*****************************************************************************
 * lxsame: return 1 if all elements of x any y match, else zero.
 *
 * Floating point versions use absolute tolerances on the allowable difference.
 *
 * $Id$
 *****************************************************************************/

#include <math.h>
#include <femdef.h>

#define EPSSP   6.0e-7
#define EPSDP   6.0e-14

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif


integer lisame (integer n,
		const integer* x, integer incx,
		const integer* y, integer incy)
{ 
  register integer i;

  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;

  for (i = 0; i < n; i++) if (x[i*incx] != y[i*incy]) return 0;
  
  return 1;
}


integer ldsame (integer n,
		const double* x, integer incx,
		const double* y, integer incy)
{ 
  register integer i;

  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;

  for (i = 0; i < n; i++) if (fabs (x[i*incx] - y[i*incy]) > EPSDP) return 0;

  return 1;
}


integer lssame (integer n,
		const float* x, integer incx,
		const float* y, integer incy)
{ 
  register integer i;

  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;

#if defined(__GNUC__) || defined(__uxp__) || defined(_SX)
  for (i = 0; i < n; i++) if (fabs  (x[i*incx] - y[i*incy]) > EPSSP) return 0;
#else
  for (i = 0; i < n; i++) if (fabsf (x[i*incx] - y[i*incy]) > EPSSP) return 0;
#endif
  
  return 1;
}
