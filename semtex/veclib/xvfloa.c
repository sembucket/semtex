/*****************************************************************************
 * Cast integer vector to floating point.                                    *
 *****************************************************************************/


void dvfloa(int n, const int *x, int incx, double *y, int incy)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = (double) *x;
    x += incx;
    y += incy;
  }
}





void svfloa(int n, const int *x, int incx, float *y, int incy)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = (float) *x;
    x += incx;
    y += incy;
  }
}
