/*****************************************************************************
 * xvlog:  y[i] = log(x[i]).
 *****************************************************************************/

#include <math.h>


void dvlog (int n, const double *x, int incx,
                         double *y, int incy)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = log (*x);
    x += incx;
    y += incy;
  }
}


void svlog (int n, const float *x, int incx,
                          float *y, int incy)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
#ifdef __GNUC__
    *y = (float) log (*x);
#else
    *y = logf (*x);
#endif
    x += incx;
    y += incy;
  }
}
