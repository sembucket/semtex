/*****************************************************************************
 * xvabs:  y[i] = abs(x[i]).
 *****************************************************************************/

#include <math.h>


void dvabs (int n, const double *x, int incx,
                         double *y, int incy)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = fabs( *x );
    x += incx;
    y += incy;
  }
}


void ivabs(int n, const int *x, int incx,
                        int *y, int incy)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = (*x<0) ? -*x : *x;
    x += incx;
    y += incy;
  }
}


void svabs (int n, const float *x, int incx,
                         float *y, int incy)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
#ifdef __GNUC__
    *y = (float) fabs (*x);
#else
    *y = fabsf (*x);
#endif
    x += incx;
    y += incy;
  }
}
