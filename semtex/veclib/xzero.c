/*****************************************************************************
 * xzero:  x[i] = 0.
 *
 * $Id$
 *****************************************************************************/

#include <string.h>
#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dzero (integer n, double* x, integer incx)
{
#if defined(DEBUG)
  register integer i;

  x += (incx < 0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = 0.0;
#else
  if (incx == 1)
    memset (x, '\0', n * sizeof (double));

  else {
    register integer i;

    x += (incx < 0) ? (-n+1)*incx : 0;

    for (i = 0; i < n; i++) x[i*incx] = 0.0;
  }
#endif
}


void izero (integer n, integer* x, integer incx)
{
  if (incx == 1)
    memset(x, '\0', n * sizeof (integer));

  else {
    register integer i;

    x += (incx < 0) ? (-n+1)*incx : 0;

    for (i = 0; i < n; i++) x[i*incx] = 0;
  }
}


void szero (integer n, float* x, integer incx)
{
  if (incx == 1)
    memset(x, '\0', n * sizeof (float));

  else {
    register integer i;

    x += (incx < 0) ? (-n+1)*incx : 0;

    for (i = 0; i < n; i++) x[i*incx] = 0.0F;
  }
}
