/*****************************************************************************
 * xvhypot:  z[i] = sqrt(SQR(x[i]) + SQR(y[i])).
 *****************************************************************************/

#include <math.h>


void dvhypot (int n, const double *x, int incx, const double *y, int incy,
	                                              double *z, int incz)
{
  register int i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) {
    *z = hypot (*x, *y);
    x += incx;
    y += incy;
    z += incz;
  }
}


void svhypot (int n, const float *x, int incx, const float *y, int incy,
	                                             float *z, int incz)
{
  register int i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) {
#if (defined(__GNUC__) || defined (__uxp__)) /* -- No single-precision maths */
    *z = (float) hypot (*x, *y);
#else
    *z = hypotf (*x, *y);
#endif
    x += incx;
    y += incy;
    z += incz;
  }
}
