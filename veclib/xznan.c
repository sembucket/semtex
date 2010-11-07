/*****************************************************************************
 * xznan:  x[i] = isnan(x[i]) ? 0.0 : x[i]
 *
 * $Id $
 *****************************************************************************/

#include <math.h>
#include <cfemdef.h>


#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dznan (integer n, double* x, integer incx)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;

#if defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) x[i*incx] = isnan (x[i*incx]) ? 0.0 : x[i*incx];
}


void sznan (integer n, float* x, integer incx)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;

#if defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) x[i*incx] = isnan (x[i*incx]) ? 0.0 : x[i*incx];
}
