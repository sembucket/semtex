/*****************************************************************************
 * xcopy:  y[i] = x[i].                                                      *
 *                                                                           *
 * NB: dcopy and scopy are also in the BLAS.                                 *
 *****************************************************************************/


void  dcopy(int n, const double *x, int incx,
                         double *y, int incy)
{
  register int  i;
  

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = *x;
    x += incx;
    y += incy;
  }
}





void  icopy(int n, const int *x, int incx,
                         int *y, int incy)
{
  register int  i;
  

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = *x;
    x += incx;
    y += incy;
  }
}





void  scopy(int n, const float *x, int incx,
                         float *y, int incy)
{
  register int  i;
  

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i=0; i<n; i++) {
    *y = *x;
    x += incx;
    y += incy;
  }
}
