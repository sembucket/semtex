/*****************************************************************************
 *  z = (w * x) - y.
 *
 * $Id$
 *****************************************************************************/

#ifdef __uxp__
#pragma global novrec
#pragma global noalias
#endif


void dvvtvm (int n, const double *w, int incw,
	            const double *x, int incx,
	            const double *y, int incy,
	                  double *z, int incz)
{
  register int i;

  w += (incw<0) ? (-n+1)*incw : 0;
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = w[i*incw] * x[i*incx] - y[i*incy];
}


void svvtvm (int n, const float *w, int incw,
	            const float *x, int incx,
	            const float *y, int incy,
	                  float *z, int incz)
{
  register int i;

  w += (incw<0) ? (-n+1)*incw : 0;
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) z[i*incz] = w[i*incw] * x[i*incx] - y[i*incy];
}
