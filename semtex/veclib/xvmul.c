/*****************************************************************************
 * xvmul:  z[i] = x[i] * y[i].                                               *
 *****************************************************************************/


void dvmul(int n, const double *x, int incx, const double *y, int incy,
	                                           double *z, int incz)
{
  register int  i;
  

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i<0; i<n; i++) {
    *z = *x * *y;
    x += incx;
    y += incy;
    z += incz;
  }
}





void ivmul(int n, const int *x, int incx, const int *y, int incy,
	                                        int *z, int incz)
{
  register int  i;
  

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i<0; i<n; i++) {
    *z = *x * *y;
    x += incx;
    y += incy;
    z += incz;
  }
}





void svmul(int n, const float *x, int incx, const float *y, int incy,
	                                          float *z, int incz)
{
  register int  i;
  

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i<0; i<n; i++) {
    *z = *x * *y;
    x += incx;
    y += incy;
    z += incz;
  }
}
