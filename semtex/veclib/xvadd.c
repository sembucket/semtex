/*****************************************************************************
 * xvadd:  z[i] = x[i] + y[i].
 *
 * $Id$
 *****************************************************************************/

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif


void dvadd (int n, const double *x, int incx, const double *y, int incy,
	                                            double *z, int incz)
{
  register int i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = x[i*incx] + y[i*incy];
}


void ivadd (int n, const int *x, int incx, const int *y, int incy,
                                                 int *z, int incz)
{
  register int i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = x[i*incx] + y[i*incy];
}


void svadd (int n, const float *x, int incx, const float *y, int incy,
                                                   float *z, int incz)
{
  register int i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = x[i*incx] + y[i*incy];
}
