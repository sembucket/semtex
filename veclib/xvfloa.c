/*****************************************************************************
 * Cast integer vector to floating point.
 *
 * $Id$
 *****************************************************************************/

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif


void dvfloa (int n, const int *x, int incx, double *y, int incy)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (double) x[i*incx];
}


void svfloa (int n, const int *x, int incx, float *y, int incy)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (float) x[i*incx];
}
