/*****************************************************************************
 * xvcos:  y[i] = cos(x[i]).
 *****************************************************************************/

#include <math.h>


void dvcos (int n, const double *x, int incx,
                         double *y, int incy)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) {
    *y = cos (*x);
    x += incx;
    y += incy;
  }
}


void svcos (int n, const float *x, int incx,
                         float *y, int incy)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) {
#ifdef __GNUC__			/* -- No support for single-precision maths. */
    *y = (float) cos (*x);
#else 
    *y = cosf (*x);
#endif
    x += incx;
    y += incy;
  }
}
