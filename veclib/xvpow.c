/*****************************************************************************
 * xvpow:  z[i] = pow(x[i], y[i]).
 *****************************************************************************/

#include <math.h>


void dvpow (int n, const double *x, int incx, const double *y, int incy,
	                                            double *z, int incz)
{
  register int     i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) {
    *z = pow (*x, *y);
    x += incx;
    y += incy;
    z += incz;
  }
}


void svpow (int n, const float *x, int incx, const float *y, int incy,
	                                           float *z, int incz)
{
  register int    i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) {
#ifdef __GNUC__
    *z = (float) pow (*x, *y);
#else
    *z = powf (*x, *y);
#endif
    x += incx;
    y += incy;
    z += incz;
  }
}
