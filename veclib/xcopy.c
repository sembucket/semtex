/*****************************************************************************
 * xcopy:  y[i] = x[i].
 *
 * NB: dcopy and scopy are also in the BLAS, but we give FORTRAN-callable
 * interface too to override these.
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


void dcopy_ (int *n, const double *x, int *incx,
                           double *y, int *incy)
{
  if ((*incx) == 1 && (*incy) == 1)
    memcpy (y, x, (*n) * sizeof (double));

  else {
    register int i;
    register int nn = *n;
    register int ix = *incx;
    register int iy = *incy;
    
    x += (ix < 0) ? (-nn + 1)*ix : 0;
    y += (iy < 0) ? (-nn + 1)*iy : 0;

    for (i = 0; i < nn; i++) y[i*iy] = x[i*ix];
  }
}


void icopy_ (int *n, const int *x, int *incx,
                           int *y, int *incy)
{
  if ((*incx) == 1 && (*incy) == 1)
    memcpy (y, x, (*n) * sizeof (int));

  else {
    register int i;
    register int nn = *n;
    register int ix = *incx;
    register int iy = *incy;
    
    x += (ix < 0) ? (-nn + 1)*ix : 0;
    y += (iy < 0) ? (-nn + 1)*iy : 0;

    for (i = 0; i < nn; i++) y[i*iy] = x[i*ix];
  }
}


void scopy_ (int *n, const float *x, int *incx,
                           float *y, int *incy)
{
  if ((*incx) == 1 && (*incy) == 1)
    memcpy (y, x, (*n) * sizeof (float));

  else {
    register int i;
    register int nn = *n;
    register int ix = *incx;
    register int iy = *incy;
    
    x += (ix < 0) ? (-nn + 1)*ix : 0;
    y += (iy < 0) ? (-nn + 1)*iy : 0;

    for (i = 0; i < nn; i++) y[i*iy] = x[i*ix];
  }
}
