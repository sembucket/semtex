/*****************************************************************************
 * xvatn2:  z[i] = atan2(x[i], y[i]).
 *****************************************************************************/

#include <math.h>


void dvatn2 (int n, const double *x, int incx, const double *y, int incy,
	                                             double *z, int incz)
{
  register int     i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) {
    *z = atan2 (*x, *y);
    x += incx;
    y += incy;
    z += incz;
  }
}


void svatn2 (int n, const float *x, int incx, const float *y, int incy,
	                                            float *z, int incz)
{
  register int    i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i = 0; i < n; i++) {
#ifdef __GNUC__
    *z = (float) atan2 (*x, *y);
#else
    *z = atan2f (*x, *y);
#endif
    x += incx;
    y += incy;
    z += incz;
  }
}
