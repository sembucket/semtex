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
  register integer  i;

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
  register integer  i;

  if (incx == 1 && incy == 1)
    memcpy (y, x, n * sizeof (float));

  else {
    x += (incx < 0) ? (-n + 1)*incx : 0;
    y += (incy < 0) ? (-n + 1)*incy : 0;

    for (i = 0; i < n; i++) y[i*incy] = x[i*incx];
  }
}

/* -- FORTRAN-callable versions to replace BLAS routines. */

void dcopy_ (integer* n, const double* x, integer* incx,
                               double* y, integer* incy)
{
  if ((*incx) == 1 && (*incy) == 1)
    memcpy (y, x, (*n) * sizeof (double));

  else {
    register integer i;
    register integer nn = *n;
    register integer ix = *incx;
    register integer iy = *incy;
    
    x += (ix < 0) ? (-nn + 1)*ix : 0;
    y += (iy < 0) ? (-nn + 1)*iy : 0;

    for (i = 0; i < nn; i++) y[i*iy] = x[i*ix];
  }
}


void icopy_ (integer* n, const integer* x, integer* incx,
                               integer* y, integer* incy)
{
  if ((*incx) == 1 && (*incy) == 1)
    memcpy (y, x, (*n) * sizeof (integer));

  else {
    register integer i;
    register integer nn = *n;
    register integer ix = *incx;
    register integer iy = *incy;
    
    x += (ix < 0) ? (-nn + 1)*ix : 0;
    y += (iy < 0) ? (-nn + 1)*iy : 0;

    for (i = 0; i < nn; i++) y[i*iy] = x[i*ix];
  }
}


void scopy_ (integer* n, const float* x, integer* incx,
                               float* y, integer* incy)
{
  if ((*incx) == 1 && (*incy) == 1)
    memcpy (y, x, (*n) * sizeof (float));

  else {
    register integer i;
    register integer nn = *n;
    register integer ix = *incx;
    register integer iy = *incy;
    
    x += (ix < 0) ? (-nn + 1)*ix : 0;
    y += (iy < 0) ? (-nn + 1)*iy : 0;

    for (i = 0; i < nn; i++) y[i*iy] = x[i*ix];
  }
}
