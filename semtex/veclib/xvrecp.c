/*****************************************************************************
 * xvrecp:  y[i] = 1.0 / x[i].
 *
 * $Id$
 *****************************************************************************/

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif


void dvrecp (int n, const double *x, int incx,
                          double *y, int incy)
{
  register int i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = 1.0 / x[i*incx];
}


void svrecp (int n, const float *x, int incx,
                          float *y, int incy)
{
  register int i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = 1.0F / x[i*incx];
}
