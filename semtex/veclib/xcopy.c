/*****************************************************************************
 * xcopy:  y[i] = x[i].
 *
 * Use memcpy for cases where both skips are unity.
 *
 * $Id$
 *****************************************************************************/

#include <string.h>
#include <femdef.h>


#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dcopy (integer n, const double* x, integer incx,
                             double* y, integer incy)
{
  register integer i;

  if (incx == 1 && incy == 1)
    memcpy (y, x, n * sizeof (double));

  else {
    x += (incx < 0) ? (-n + 1)*incx : 0;
    y += (incy < 0) ? (-n + 1)*incy : 0;

    for (i = 0; i < n; i++) y[i*incy] = x[i*incx];
  }
}


void icopy (integer n, const integer* x, integer incx,
                             integer* y, integer incy)
{
  register integer i;

  if (incx == 1 && incy == 1)
    memcpy (y, x, n * sizeof (integer));

  else {
    x += (incx < 0) ? (-n + 1)*incx : 0;
    y += (incy < 0) ? (-n + 1)*incy : 0;

    for (i = 0; i < n; i++) y[i*incy] = x[i*incx];
  }
}


void scopy (integer n, const float* x, integer incx,
                             float* y, integer incy)
{
  register integer i;

  if (incx == 1 && incy == 1)
    memcpy (y, x, n * sizeof (float));

  else {
    x += (incx < 0) ? (-n + 1)*incx : 0;
    y += (incy < 0) ? (-n + 1)*incy : 0;

    for (i = 0; i < n; i++) y[i*incy] = x[i*incx];
  }
}
