/*****************************************************************************
 * xvexp:  y[i] = exp (x[i]).
 *****************************************************************************/

#include <math.h>


void dvexp (int n, const double *x, int incx,
                         double *y, int incy)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = exp(*x);
    x += incx;
    y += incy;
  }
}


void svexp (int n, const float *x, int incx,
                         float *y, int incy)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
#ifdef __GNUC__
    *y = (float) exp (*x);
#else
    *y = expf (*x);
#endif
    x += incx;
    y += incy;
  }
}





