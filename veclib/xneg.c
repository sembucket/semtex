/*****************************************************************************
 * xneg:  x[i] = -x[i].
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>


void dneg (integer n, double* x, integer incx)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;

#if defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) x[i*incx] = -x[i*incx];
}


void ineg (integer n, integer* x, integer incx)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;

#if defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) x[i*incx] = -x[i*incx];
}


void sneg (integer n, float* x, integer incx)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;

#if defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) x[i*incx] = -x[i*incx];
}
