/*****************************************************************************
 * xmask:  conditional assignment:  if (y[i]) z[i] = w[i]; else z[i] = x[i].
 *
 * $Id$
 *****************************************************************************/


void dmask (int n, const double* w, int incw, const double* x, int incx,
	           const int*    y, int incy,       double* z, int incz)
{
  register int  i;

  w += (incw < 0) ? (-n + 1) * incw : 0;
  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;
  z += (incz < 0) ? (-n + 1) * incz : 0;

  for (i = 0; i < n; i++) {
    *z = (*y) ? *w : *x;
    w += incw;
    x += incx;
    y += incy;
    z += incz;
  }
}


void imask (int n, const int* w, int incw, const int* x, int incx,
	           const int* y, int incy,       int* z, int incz)
{
  register int  i;

  w += (incw < 0) ? (-n + 1) * incw : 0;
  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;
  z += (incz < 0) ? (-n + 1) * incz : 0;

  for (i = 0; i < n; i++) {
    *z = (*y) ? *w : *x;
    w += incw;
    x += incx;
    y += incy;
    z += incz;
  }
}


void smask (int n, const float* w, int incw, const float* x, int incx,
	           const int*   y, int incy,       float* z, int incz)
{
  register int  i;

  w += (incw < 0) ? (-n + 1) * incw : 0;
  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;
  z += (incz < 0) ? (-n + 1) * incz : 0;

  for (i = 0; i < n; i++) {
    *z = (*y) ? *w : *x;
    w += incw;
    x += incx;
    y += incy;
    z += incz;
  }
}
