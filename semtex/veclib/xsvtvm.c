/*****************************************************************************
 * xsvtvm:   z[i] = alpha * x[i] - y[i].                                     *
 *****************************************************************************/


void dsvtvm(int n, double alpha, const double *x, int incx,
                                 const double *y, int incy,
                                       double *z, int incz)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i=0; i<n; i++) {
    *z = alpha * *x - *y;
    x += incx;
    y += incy;
    z += incz;
  }
}





void  ssvtvm(int n, float alpha, const float *x, int incx,
                                 const float *y, int incy,
                                       float *z, int incz)
{
  register int  i;


  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;
  z += (incz<0) ? (-n+1)*incz : 0;

  for (i=0; i<n; i++) {
    *z = alpha * *x - *y;
    x += incx;
    y += incy;
    z += incz;
  }
}
