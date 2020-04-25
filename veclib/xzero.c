/*****************************************************************************
 * xzero:  x[i] = 0.
 *
 * $Id: xzero.c,v 9.1 2019/05/30 06:36:14 hmb Exp $
 *****************************************************************************/

#include <string.h>
#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dzero (int_t n, double* x, int_t incx)
{
#if defined(DEBUG)
  register int_t i;

  x += (incx < 0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = 0.0;
#else
  if (incx == 1)
    memset (x, '\0', n * sizeof (double));

  else {
    register int_t i;

    x += (incx < 0) ? (-n+1)*incx : 0;

    for (i = 0; i < n; i++) x[i*incx] = 0.0;
  }
#endif
}


void izero (int_t n, int_t* x, int_t incx)
{
  if (incx == 1)
    memset(x, '\0', n * sizeof (int_t));

  else {
    register int_t i;

    x += (incx < 0) ? (-n+1)*incx : 0;

    for (i = 0; i < n; i++) x[i*incx] = 0;
  }
}


void szero (int_t n, float* x, int_t incx)
{
  if (incx == 1)
    memset(x, '\0', n * sizeof (float));

  else {
    register int_t i;

    x += (incx < 0) ? (-n+1)*incx : 0;

    for (i = 0; i < n; i++) x[i*incx] = 0.0F;
  }
}
