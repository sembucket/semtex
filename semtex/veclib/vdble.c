void vdble(int n, const float *x, int incx, double *y, int incy)
/* ========================================================================= *
 * Cast single to double precision.                                          *
 * ========================================================================= */
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
