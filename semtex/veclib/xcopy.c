/*****************************************************************************
 * xcopy:  y[i] = x[i].
 *
 * NB: dcopy and scopy are also in the BLAS.
 * We do fast copy for cases where both skips are unity.
 *
 * $Id$
 *****************************************************************************/

#include <string.h>

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif


void dcopy (int n, const double *x, int incx,
                         double *y, int incy)
{
  register int i;

  if (incx == 1 && incy == 1)
    memcpy (y, x, n * sizeof (double));

  else {
    x += (incx < 0) ? (-n + 1)*incx : 0;
    y += (incy < 0) ? (-n + 1)*incy : 0;

    for (i = 0; i < n; i++) y[i*incy] = x[i*incx];
  }
}


void icopy (int n, const int *x, int incx,
                         int *y, int incy)
{
  register int  i;

  if (incx == 1 && incy == 1)
    memcpy (y, x, n * sizeof (int));

  else {
    x += (incx < 0) ? (-n + 1)*incx : 0;
    y += (incy < 0) ? (-n + 1)*incy : 0;

    for (i = 0; i < n; i++) y[i*incy] = x[i*incx];
  }
}





void scopy (int n, const float *x, int incx,
                         float *y, int incy)
{
  register int  i;

  if (incx == 1 && incy == 1)
    memcpy (y, x, n * sizeof (float));

  else {
    x += (incx < 0) ? (-n + 1)*incx : 0;
    y += (incy < 0) ? (-n + 1)*incy : 0;

    for (i = 0; i < n; i++) y[i*incy] = x[i*incx];
  }
}
