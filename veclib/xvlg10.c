/*****************************************************************************
 * xvlg10:  y[i] = log10(x[i]).
 *****************************************************************************/

#include <math.h>


void dvlg10 (int n, const double *x, int incx,
                          double *y, int incy)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = log10 (*x);
    x += incx;
    y += incy;
  }
}


void svlg10 (int n, const float *x, int incx,
                          float *y, int incy)
{
  register int  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
#ifdef __GNUC__
    *y = (float) log10 (*x);
#else
    *y = log10f (*x);
#endif
    x += incx;
    y += incy;
  }
}
