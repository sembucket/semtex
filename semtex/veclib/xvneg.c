/*****************************************************************************
 * xvneg:  y[i] = -x[i].
 *
 * $Id$
 *****************************************************************************/

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif


void dvneg (int n, const double *x, int incx,
                         double *y, int incy)
{
  register int i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = -x[i*incx];
}


void ivneg (int n, const int *x, int incx,
                         int *y, int incy)
{
  register int i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = -x[i*incx];
}


void svneg (int n, const float *x, int incx,
                         float *y, int incy)
{
  register int i;
  
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = -x[i*incx];
}
