/*****************************************************************************
 * xvneg:  y[i] = -x[i].                                                     *
 *****************************************************************************/


void  dvneg(int n, const double *x, int incx,
                         double *y, int incy)
{
  register int  i;
  

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = -*x;
    x += incx;
    y += incy;
  }
}





void  ivneg(int n, const int *x, int incx,
                         int *y, int incy)
{
  register int  i;
  

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = -*x;
    x += incx;
    y += incy;
  }
}





void  svneg(int n, const float *x, int incx,
                         float *y, int incy)
{
  register int  i;
  

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = -*x;
    x += incx;
    y += incy;
  }
}
