/*****************************************************************************
 * Cast double to single precision.                                          *
 *****************************************************************************/


void vsngl(int n, const double *x, int incx, float *y, int incy)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y  = (float) *x;
    x += incx;
    y += incy;
  }
}
