/*****************************************************************************
 * xvsqrt:  y[i] = sqrt(x[i]).
 *****************************************************************************/

#include <math.h>


void dvsqrt (int n, const double *x, int incx,
                          double *y, int incy)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = sqrt (*x);
    x += incx;
    y += incy;
  }
}


void svsqrt (int n, const float *x, int incx,
                          float *y, int incy)
{
  register int  i;
  
  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
#ifdef __GNUC__			/* -- No support for single-precision maths. */
    *y = (float) sqrt (*x);
#else
    *y = sqrtf (*x);
#endif
    x += incx;
    y += incy;
  }
}
