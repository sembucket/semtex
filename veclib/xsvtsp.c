/*****************************************************************************
 * xsvtsp:  y[i] = alpha * x[i] + beta.
 *
 * $Id$
 *****************************************************************************/

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif


void dsvtsp (int n, double alpha, double beta, const double *x, int incx,
                                                     double *y, int incy)
{
  register int  i;
  
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = alpha * x[i*incx] + beta;
}


void svtsp (int n, float alpha, float beta, const float *x, int incx,
                                                  float *y, int incy)
{
  register int  i;
  
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = alpha * x[i*incx] + beta;
}
