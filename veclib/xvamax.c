/*****************************************************************************
 * xvamax:  z[i] = MAX(ABS(x[i]), ABS(y[i])).
 *****************************************************************************/

#include <stdlib.h>
#include <math.h>

#define MAX(x, y) ( ((x)>(y)) ? (x) : (y))


void dvamax (int n, const double *x, int incx, const double *y, int incy,
	                                             double *z, int incz)
{
  register int     i;
  register double  absx, absy;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i=0; i<n; i++) {
    absx = fabs (*x);
    absy = fabs (*y);
    *z = MAX (absx, absy);
    x += incx;
    y += incy;
    z += incz;
  }
}


void ivamax (int n, const int *x, int incx, const int *y, int incy,
                                                  int *z, int incz)
{
  register int i, absx, absy;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i=0; i<n; i++) {
    absx = abs (*x);
    absy = abs (*y);
    *z = MAX (absx, absy);
    x += incx;
    y += incy;
    z += incz;
  }
}


void svamax (int n, const float *x, int incx, const float *y, int incy,
	                                            float *z, int incz)
{
  register int    i;
  register float  absx, absy;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i=0; i<n; i++) {
#ifdef __GNUC__
    absx = (float) fabs (*x);
    absy = (float) fabs (*y);
#else
    absx = fabsf (*x);
    absy = fabsf (*y);
#endif
    *z = MAX (absx, absy);
    x += incx;
    y += incy;
    z += incz;
  }
}
