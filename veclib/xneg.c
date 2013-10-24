/*****************************************************************************
 * xneg:  x[i] = -x[i].
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef.h>


#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dneg (int_t n, double* x, int_t incx)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;

#if defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) x[i*incx] = -x[i*incx];
}


void ineg (int_t n, int_t* x, int_t incx)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;

#if defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) x[i*incx] = -x[i*incx];
}


void sneg (int_t n, float* x, int_t incx)
{
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;

#if defined(_SX)
#pragma vdir nodep
#endif

  for (i = 0; i < n; i++) x[i*incx] = -x[i*incx];
}
