/*****************************************************************************
 * xssub:  y[i] = alpha - x[i].                                              *
 *****************************************************************************/


void dssub(int n, double alpha, const double *x, int incx,
	                              double *y, int incy)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = alpha - *x;
    x += incx;
    y += incy;
  }
}





void isub(int n, int alpha, const int *x, int incx,
	                          int *y, int incy)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = alpha - *x;
    x += incx;
    y += incy;
  }
}





void sssub(int n, float alpha, const float *x, int incx,
	                             float *y, int incy)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = alpha - *x;
    x += incx;
    y += incy;
  }
}
